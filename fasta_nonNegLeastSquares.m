%   Solve the non-negative least-squares problem
%           min ||Ax-b||^2
%           subject to x>=0
%   using the solver FASTA.  
%
%  Inputs:
%    A   : A matrix or function handle
%    At  : The adjoint/transpose of A
%    b   : A column vector of measurements
%    x0  : Initial guess of solution, often just a vector of zeros
%    opts: Optional inputs to FASTA
%
%   For this code to run, the solver "fasta.m" must be in your path.
%
%   For more details, see the FASTA user guide, or the paper "A field guide
%   to forward-backward splitting with a FASTA implementation."
%
%   Copyright: Tom Goldstein, 2014.



function [ solution, outs ] = fasta_nonNegLeastSquares( A,At,b,mu,x0,opts,l1 )

%%  Check whether we have function handles or matrices
if ~isnumeric(A)
    assert(~isnumeric(At),'If A is a function handle, then At must be a handle as well.')
end
%  If we have matrices, create handles just to keep things uniform below
if isnumeric(A)
    At = @(x)A'*x;
    A = @(x) A*x;
end

%  Check for 'opts'  struct
if ~exist('opts','var') % if user didn't pass this arg, then create it
    opts = [];
end


%%  Define ingredients for FASTA
%  Note: fasta solves min f(Ax)+g(x).
%  f(z) = .5 ||z - b||^2
f    = @(z) .5*norm(z-b,'fro')^2;
grad = @(z) z-b;
% g(z) = mu*|z|
g = @(x) 0;
% proxg(z,t) = argmin .5||x-z||^2, subject to x>=0
if l1
    prox = @(x,t) max(x,0) + project_1ball(x,mu);
else
     prox = @(x,t) max(x,0);
end

%% Call solver
[solution, outs] = fasta(A,At,f,grad,g,prox,x0,opts);

end
function [ z ] = project_1ball( z,radius )
%  By Moreau's identity, projection onto 1-norm ball can be computed 
%  using the proximal of the conjugate problem, which is L-infinity
%  minimization.
  z = z -  prox_infinityNorm(z,radius);    
end


%% Perform prox operator:   min ||x||_inf + (1/2t)||x-w||^2
function [ xk ] = prox_infinityNorm( w,t )
    N = length(w);
    wabs = abs(w);
    ws = (cumsum(sort(wabs,'descend'))- t)./(1:N)';
    alphaopt = max(ws);
    if alphaopt>0 
      xk = min(wabs,alphaopt).*sign(w); % truncation step
    else
      xk = zeros(size(w)); % if t is big, then solution is zero
    end       
end