close all;
clear all;

data_folder = '../../../datasets/';
load([data_folder 'netnmfsc_comparison.mat'])
load('resampled_netnmfsc.mat')

%%
ix=1;
r=5;

% set stable options
opts.LQF = 10;
opts.l1 = 0;
opts.randInit=false;
opts.fasta.recordObjective=false;
opts.smoothInit=true;
opts.fasta.tol=1e-10;
opts.fasta.maxIters = 1000;
opts.maxIters=500;
opts.earlyStop=false;
opts.fasta.verbose=false;
opts.fasta.adaptive=true;
opts.fasta.accelerate=true;
opts.fasta.tau = 10;
opts.netnmfscStop=false;
opts.verbose=true;
z = library_size_normalization(X1{ix});
[Y,W,H,obj,W0,H0] = boxR2RNNGLS(z,A,r,opts);
%%

save('smooth_initializations', 'W0','H0')
%%

L = diag(sum(A,1)) - A;
mask = (z==0);
WH = W0*H0';
norm(z(mask)-WH(mask),'F')^2 + trace(W0'*L*W0)