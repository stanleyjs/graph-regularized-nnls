%%  This script evaluates FASTA performance under different solver options
% 1) Load example data from netnmfsc
% 2) Create bigger datasets by multinomial resampling
% 3) Test time, mse, etc for fasta defaults, no adaptation, acceleration,
% and acceleration + adaptation.
% Feb 8 2021: RandIndex appears to have some bugs in processing.
close all;
clear all;
load('netnmfsc_comparison.mat')
X = double(X);
[m,n] = size(X);
netnmfsc.labels = double(clusters)+1;
netnmfsc.Y = double(Wpy*Hpy);
netnmfsc.genes = netgenes;
netnmfsc.H = double(Hpy);
netnmfsc.W = double(Wpy);

r = 5;
%% 2) add some more cells for time-based comparison. Simple model is to resample the clusters.
close all;
clear all;
load('netnmfsc_comparison.mat')

median_library_size = ceil(median(sum(X,1)));
k = length(unique(netnmfsc.labels));
cluster_size = zeros(1,k);
pdist = zeros(m,k);
X1{1} = [];
labels1{1} = [];
for cluster = 1:k
    cells = netnmfsc.labels==cluster;
    cluster_size = nnz(cells);
    pdist(:,cluster) = sum(X(:,cells),2);
    pdist(:,cluster) = pdist(:,cluster)/sum(pdist(:,cluster));
    X1{1} = [X1{1} mnrnd(median_library_size,pdist(:,cluster),cluster_size)'];
    labels1{1} = [labels1{1}; cluster*ones(cluster_size,1);];
end

for multiple= 2:5
    X1{multiple}=[];
    labels1{multiple}=[];
    for cluster=1:k
        cells = labels1{multiple-1}==cluster;
        cluster_size = nnz(cells);
        new_member_size = cluster_size;
        new_members = mnrnd(median_library_size,pdist(:,cluster),new_member_size)';
        X1{multiple} = [X1{multiple} X1{multiple-1}(:,cells) ...
            new_members];
        labels1{multiple} = [labels1{multiple}; cluster*ones(cluster_size+new_member_size,1)];
    end
end





%% 3) run experiment
clear results

% tables for results
experiments = {'defaults', 'noadapt','accelerate+noadapt','adaptive+accelerate'};
results = struct();
results.mse = {[],[],[],[]};
results.gradientnorm ={[],[],[],[]};
results.obj = {[],[],[],[]};
results.time = {[],[],[],[]};
results.silhouette = {[],[],[],[]};
results.ari = {[],[],[],[]};
% experiments
for multiple = 1:5
    disp(['current multiple ', num2str(multiple)])
    clear opts %reset to defaults
    
    % set stable options
    opts.LQF = 10;
    opts.l1 = 0;
    opts.randInit=false;
    opts.fasta.recordObjective=false;
    opts.randInit=false;
    opts.smoothInit=true;
    opts.fasta.tol=1e-6;
    opts.maxIters=40;
    opts.fasta.verbose=true;
    %construct initializations
    z = library_size_normalization(X1{multiple});
    
    [opts.initialization.W,opts.initialization.H] = nnmf(z,r);
    for id = 1:numel(experiments)
       experiment = experiments{id};
       switch experiment
           case 'defaults'
               
           case 'noadapt'
               opts.fasta.adaptive = false;
               opts.fasta.accelerate = false;
           case 'accelerate+noadapt'
               opts.fasta.adaptive = false;
               opts.fasta.accelerate = true;
           case 'adaptive+accelerate'
               opts.fasta.adaptive = true;
               opts.fasta.accelerate = true;
       end
       
       
       [Y,W,H,obj] = boxR2RNNGLS(z,A,r,opts);
       
       %record time, rmse, objective results
       t_update = [obj.time.init obj.time.iter];
       sdif = opts.maxIters - size(t_update,2);
       if sdif>0
           pad = zeros(1,sdif);
           t_update = [t_update,pad];
           obj.mse = [obj.mse,pad];
           obj.obj = [obj.obj, pad];
           obj.gradientnorm = [obj.gradientnorm, pad];
       end
       results.time{id} = [results.time{id} t_update'];
       results.mse{id} = [results.mse{id} obj.mse'];
       results.obj{id} = [results.obj{id} obj.obj'];
       results.gradientnorm{id} = [results.gradientnorm{id} obj.gradientnorm'];
       %compute clustering results
       [c,labs,score] = select_clusters(H);
       score = score(2);
       clusters = labs(:,2);
       ari = RandIndex(clusters,labels1{multiple}); %Feb 8 2021: rand index looks bugged
       %record clustering results
       results.silhouette{id} = [results.silhouette{id} score];
       results.ari{id} = [results.ari{id} ari];

    end

end


%%
save("fasta_opts_comparison")


%%

data_sizes = cellfun(@(z) size(z,2),X1,'UniformOutput', false);
data_sizes = cell2mat(data_sizes);
figure
suptitle("FASTA step options comparison")
hold on;
title("40 iters time")
for id = 1:numel(experiments)
    plot((data_sizes),log(sum(results.time{id},1)),'LineWidth',2)
end
xlabel('data size')
ylabel('log time (s)')