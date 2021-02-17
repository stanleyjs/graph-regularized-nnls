close all;
clear all;
data_folder = '../../../data/';
load([data_folder 'netnmfsc_comparison.mat'])
load('resampled_netnmfsc.mat')

%%
ix=5;
r=5;

% set stable options
opts.LQF = 10;
opts.l1 = 0;
opts.randInit=true;
opts.fasta.recordObjective=false;
opts.smoothInit=true;
opts.fasta.tol=1e-6;
opts.maxIters=2;
opts.fasta.verbose=true;
opts.fasta.adaptive=true;
opts.fasta.accelerate=false;
z = library_size_normalization(X1{ix});
[~] = boxR2RNNGLS(z,A,r,opts);
