close all;
clear all;
load('/home/jay/ExtraDrive1/External/experiments/r2rilsexperiment/laplacian_regularization/boxR2RNNGLS/test/netnmfsc_comparison.mat')
X = double(X);
[m,n] = size(X);
netnmfsc.labels = double(clusters)+1;
netnmfsc.Y = double(Wpy*Hpy);
netnmfsc.genes = netgenes;
netnmfsc.H = double(Hpy);
netnmfsc.W = double(Wpy);
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


save('resampled_netnmfsc','X1')