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
%% add some more cells for time-based comparison. Simple model is to resample the clusters.
median_library_size = ceil(median(sum(X,1)));
k = length(unique(netnmfsc.labels));
cluster_size = zeros(1,k);
pdist = zeros(m,k);
X1{1} = X;
labels1{1} = netnmfsc.labels;
for cluster = 1:k
    cells = netnmfsc.labels==cluster;
    cluster_size = nnz(cells);
    pdist(:,cluster) = sum(X(:,cells),2);
    pdist(:,cluster) = pdist(:,cluster)/sum(pdist(:,cluster));
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


%%
z = library_size_normalization(X1{1});
opts.LQF = 10;
opts.l1 = 0.001;
opts.randInit=false;
opts.fasta.recordObjective=true;
opts.fasta.tol=1e-6;
opts.fasta.verbose=false;
opts.fasta.accelerate=true;
opts.fasta.adaptive=true;
opts.fasta.restart=true;
opts.maxIters=40;
[Y,W,H,obj] = boxR2RNNGLS(z,A,r,opts);
%%
[k,labs,score] = select_clusters(netnmfsc.H');
score
embedding = tsne(netnmfsc.H');
scatter(embedding(:,1),embedding(:,2),20,labs(:,k-1),'filled');
