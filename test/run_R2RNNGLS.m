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


z = library_size_normalization(X);
opts.LQF = 10;
opts.l1 = 0;
[opts.initialization.W,opts.initialization.H] = nnmf(z,r);
opts.randInit=false;
opts.smoothInit=true;
opts.fasta.recordObjective=false;
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
