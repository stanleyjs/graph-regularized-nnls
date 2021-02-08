function [k,labs,score] = select_clusters(X)
%SELECT_CLUSTERS Summary of this function goes here
%   Detailed explanation goes here
max_clusters = 20;
cluster_range = [ 2: max_clusters];
avgs = zeros(length(cluster_range),1);
clusters = zeros(size(X,1),length(cluster_range));
for ix=2:length(cluster_range)+1
    clusters(:,ix-1) = kmeans(X, cluster_range(ix-1));
    silhouette_avg = mean(silhouette(X, clusters(:,ix-1)));
    avgs(ix-1) = (silhouette_avg);
end
[~,k] = max(avgs);
labs = clusters;
score = avgs;
% score = avgs(k);
% labs = clusters(:,k);
k = cluster_range(k);
end

