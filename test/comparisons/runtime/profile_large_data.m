close all;
clear all;
load('/home/jay/ExtraDrive1/External/experiments/r2rilsexperiment/laplacian_regularization/boxR2RNNGLS/test/netnmfsc_comparison.mat')
load('resampled_netnmfsc.mat')

%%
ix=5;
r=5;

% set stable options
opts.LQF = 10;
opts.l1 = 0;
opts.randInit=true;
opts.fasta.recordObjective=false;
opts.smoothInit=false;
opts.fasta.tol=1e-6;
opts.maxIters=2;
opts.fasta.verbose=true;
opts.fasta.adaptive=true;
opts.fasta.accelerate=false;
z = library_size_normalization(X1{ix});
[~] = boxR2RNNGLS(z,A,r,opts);
