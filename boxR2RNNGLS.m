function [Y,W,H,stats,W0,H0] = boxR2RNNGLS(X, A_w, r, opts)

    %  Inputs:
    %    X   : m x n data matrix to be factored
    %    A_w : m x m Adjacency matrix for features of X 
    %          (see 'opts.L' for passing precomputed laplacian)
    %    r   : rank of factorization
    %    opts: Optional inputs to boxR2RNNGLS
    %       opts.fasta: Optional inputs to FASTA
    
    
    if ~exist('opts','var') % if user didn't pass this arg, then create it
        opts = struct();

    end
    if ~exist('r','var')
        r = [];
    end
    if nargout > 2
        stats = struct();
        % for housekeeping, go ahead and build the fields:
        stats.mse = [];
        stats.time.init = 0;
        stats.time.iter = [];
        stats.time.final = 0;
        stats.obj = [];
        stats.gradientnorm = [];
        tic
    end

    [opts,X,A_w,D_w,L_w,r,sz] = processInputs(opts, X, A_w,r);

    m = sz(1);
    n = sz(2);
    p = (m + n) * r; % the number of variables in the least squares problem.

    %% find nonzero entries for masking the data.
    if opts.completion
        [omega(:,1),omega(:,2)] = find(X);
        nv = size(omega,1);
        x = zeros(nv,1);
        for counter = 1:nv
            x(counter) = X(omega(counter,1),omega(counter,2));
        end
    else 
        error('denoising (no mask) is not yet supported');
    end
    
    %% objectives & updates: we are using a lot of extra variables here to 
    % maybe reduce some flops from repetitive slicing by using these
    % functions carefully.
    msef = @(vA,u) 1/(nv) .*sum((x-vA*u).^2);
    du = @(z) z(1:m*r);
    dv = @(z) z(m*r+1:end);
    
%     reconstruction = @(z,A,vA,u0) 
    
    l1 = @(v) sum(abs(v));
    l2 = @(df) sum(df.^2);
    lqf = @(u) compute_LQF(u,L_w,m,r);
    
    
    obj = @(z,A,vA,f0,u0,du,dv)  reconstruction(x,A,u0,du,dv);
    obj_g  = @(z) 0;
    %gradients
    gradf = @(z,A,vA,f0,u0,du,dv) reconstruction_gradient(x,A,u0,du,dv);
    lqf_grad = @(u) [compute_LQF_gradient(u,L_w,m,r);zeros(n*r,1)];
    l2_grad = @(z) z;
    
    % prox objective in the non-l1 case.
    proxg = @(z,f0,u0,du,v0,dv) box_prox(z,f0);
    
    if opts.LQF>0
        obj = @(z,A,vA,f0,u0,du,dv) obj(z,A,vA,f0,u0,du,dv) + opts.LQF * lqf(du+u0);
        gradf = @(z,A,vA,f0,u0,du,dv) gradf(z,A,vA,f0,u0,du,dv) + opts.LQF * lqf_grad(du+u0);
    end
    if opts.l2>0
        obj = @(z,A,vA,f0,u0,du,dv) obj(z,A,vA,f0,u0,du,dv) + opts.l2 * l2(z);
        gradf = @(z,A,vA,f0,u0,du,dv) gradf(z,A,vA,f0,u0,du,dv) + opts.l2 *l2_grad(z);
    end
    if opts.l1>0
        %this prox will be totally different than the box - it is not a
        %sum of proxes.
        proxg = @(z,f0,u0,du,v0,dv)  l1_box_prox(z,f0,opts.l1,m,r); 
        obj_g = @(v) opts.l1*l1(v);
        obj = @(z,A,vA,f0,u0,du) obj(z,A,vA,f0,u0,du);
    end
    
    
        
    %% initialization
    if isstruct(opts.initialization)
        Wx = opts.initialization.W;
        Hx = opts.initialization.H;
    else
        if opts.randInit
            Wx = abs(randn(m,r));
            Hx = abs(randn(r,n));
        else
            [Wx,Hx] = nnmf(X,r);

        end
    end
    if opts.smoothInit
        [Wx,~,~] = svds(L_w,r);
        Wx = Wx-min(Wx);
        
        wx = reshape(Wx,[],1);    
        A = zeros(nv,n*r);
        %MATRIX FORM
        for counter = 1:nv
            j = omega(counter,1); k = omega(counter,2);
%             Hix = j:m:(m*r);
%             A(counter,Hix) = H(k,:);
            Wix = ((k-1)*r)+1;
            Wix2 = (Wix)+r-1;
            A(counter,Wix:Wix2) = Wx(j,:);
        end

        A = sparse(A);
        Hx = fasta_nonNegLeastSquares( A,A',x,0,zeros(n*r,1),opts.fasta,0 );
        Hx = reshape(Hx,r,n)';
    end
    
    if all(size(Hx) == [r,n])
        Hx = Hx';
    end
    
    
    if nargout>3 %% RETURN THE INITIAL MATRICES: FOR TESTING
        W0 = Wx;
        H0 = Hx;
    end
    wx = reshape(Wx,[],1);
    hx = reshape(Hx',[],1);
    maxIters = opts.maxIters;
    f = [wx;hx];
    
    fasta_A = @(x) x;
    df = zeros(p,1);
    [A0,vA0] = compute_A(f,omega,m,r,nv,p);
    obj0 = obj(zeros(p,1),A0,vA0,f,wx,zeros(m*r,1),zeros(n*r,1));
    if opts.verbose
        fprintf(['************** R2R Initialized ************** \n'])
        fprintf(['******** Objective :',num2str(obj0) ' ********\n'])
    end
    if nargout >2
        stats.time.init = toc;
        stats.obj = obj0;
        
    end
    
    %% outer loop
    for t = 1:maxIters-1
        if mod(t+1,5) ==0 & opts.verbose
            fprintf(['************** R2R Iterate:', num2str(t+1), ' ************** \n'])
            fprintf(['******** Objective :',num2str(objt(zeros(p,1))) ' ********\n'])
        end
        if nargout > 2
            tic
        end
        f = df+f;
        [A,vA] = compute_A(f,omega,m,r,nv,p);

        
        %ideally these sliced calls should be all done at once in FASTA to
        %save a few flops..
        u0 = du(f);
        v0 = dv(f);

        objt = @(z) obj(z,A,vA,f,u0,du(z),dv(z)); 
        gradt = @(z) gradf(z,A,vA,f,u0,du(z),dv(z));
        proxt = @(z,foo) proxg(z,f,u0,du(z),v0,dv(z));
        gt = @(z) obj_g(v0+dv(z));
        [df,~] = fasta(fasta_A,fasta_A,objt,gradt,gt,proxt,df,opts.fasta);
        
        if opts.earlyStop
            % these values should probably be computed from fasta
            % also consider location... this is pre update to f, so it is
            % the t-1 iterate.
            if t>1
                mse_eval0 = mse_eval;
                obj_eval0 = obj_eval;
            end
            mse_eval = msef(vA,u0);
            obj_eval = obj(zeros(p,1),A,vA,f,u0,zeros(m*r,1),zeros(n*r,1)) + ...
                obj_g(v0);
            gradient_norm = sqrt(sum(df.^2));
            
        else
            mse_eval = [];
            obj_eval = [];
            gradient_norm = [];
        end
        if nargout > 2
            if isempty(mse_eval)
                mse_eval = msef(vA,u0);
                obj_eval = obj(zeros(p,1),A,vA,f,u0,zeros(m*r,1),zeros(n*r,1)) + ...
                obj_g(v0);
                gradient_norm = sqrt(sum(df.^2));
            end
            stats.time.iter = [stats.time.iter toc];
            stats.mse = [stats.mse  msef(vA,u0)];
            % zeros are used to compute objective here because df is the
            % objective at the next step..
            stats.obj = [stats.obj obj_eval];
            stats.gradientnorm = [stats.gradientnorm gradient_norm];
            
        end
        if opts.earlyStop
            if gradient_norm/sqrt(p) <opts.tol
                %updates are too small: this is equivalent to checking 
                % the mse between the previous iterate 
                break;
            end
            if t>1
                if abs(mse_eval0/mse_eval)<opts.tol
                        % ANOTHER EARLY STOPPING: IF OBSERVED ERROR DID NOT CHANGE OVER THE
                         % ITERATIONS

                    break;
                end
                if abs(obj_eval0/obj_eval)<opts.tol
                    % objective did not change much.
                    break;
                end
            end
        end
        if opts.netnmfscStop
            if t>1
                if abs(obj_eval0-obj_eval)<opts.netnmfscStop
                    break;
                end
            end
        end
            


    end
f = df +f;
if nargout >2
    stats.mse = [stats.mse msef(vA,u0)];
    stats.obj = [stats.obj obj(zeros(p,1),A,vA,f,du(f),zeros(m*r,1),zeros(n*r,1)) + obj_g(dv(f))];

    tic
end
W = reshape(f(1:m*r),m,r);
H = reshape(f(m*r+1:end),r,n)';
Y = W*H';
if nargout > 2
    stats.time.final = toc + sum(stats.time.iter);
end
end 
function [penalty] = reconstruction(x,A,u0,du,dv)

    penalty = [du+u0;dv];
    penalty = A*penalty;
    penalty = x-penalty;
    penalty = (penalty.^2);
    penalty = sum(penalty);
end
function [gradient] = reconstruction_gradient(x,A,u0,du,dv)

    gradient = [du+u0;dv];
    gradient = A*gradient;
    gradient = x-gradient;
    gradient = -2.*A'*gradient;
    
end
function [A,vA] = compute_A(f,omega,m,r,nv,p)
    %CONSTRUCTION OF MATRIX A, 
    
    A = zeros(nv,p);
        
        %MATRIX FORM
%     for counter = 1:nv
%         j = omega(counter,1); k = omega(counter,2);
%          Hix = j:m:(m*r);
%          A(counter,Hix) = H(k,:);
%         Wix = (m*r)+((k-1)*r)+1;
%         Wix2 = (Wix)+r-1;
%         A(counter,Wix:Wix2) = W(j,:);
%     end

    for counter = 1:nv
        j = omega(counter,1); k = omega(counter,2);
        Hix = j:m:(m*r);
        A(counter,Hix) = f(m*r+(k-1)*r+1:m*r+k*r); %the k-th column of H
        Wix = (m*r)+((k-1)*r)+1;
        Wix2 = (Wix)+r-1;
        A(counter,Wix:Wix2) = f(j:m:m*r);
    end
    A = sparse(A);
    vA = sparse(A(:,1:m*r));
        
end

function [z] =l1_box_prox(z, f0, mu, m, r)
    disp("l1 shrink")
    w = box_prox(z(1:m*r),f0(1:m*r));
    h0 = f0(m*r+1:end);
    dh = z(m*r+1:end);
    h = dh +h0;
    %collect the indices to modify first before changing them to prevent double
    %modification
    pos_inds = h > mu;
    neg_inds = h < -1 * mu; %this should never be reached
    zero_inds = abs(h) <= mu;

    h(pos_inds) = max(h(pos_inds)-mu,0) - h0(pos_inds);
    h(neg_inds) = max(h(neg_inds)+mu,0) -  h0(neg_inds);
    h(zero_inds) = 0 - h0(zero_inds);

    z =  [w;h];
end


function [z] = box_prox(z,f0)
    z((z+f0)<0) = 0;
    
end

function penalty = compute_LQF(z,L,m,r)

%try the big kronecker product: feb 7 2020
% requires changing L in the initialization to be L = kron(speye(r),L);

% penalty = z((1:m*r))'*L;
% penalty = penalty * z(1:m*r);

%%%
%the following appears to be faster than just multiplying by the big kronecker product.
penalty = 0;
    for j = 1:r
        Lv = L*z((j-1)*m+1:j*m);
        penalty = penalty + z((j-1)*m+1:j*m)'*Lv;
    end
end

function update = compute_LQF_gradient(z,L,m,r)  

%try the big kronecker product: feb 7 2020
% requires changing L in the initialization to be L = kron(speye(r),L);
% update = L*z((1:m*r));

%%%
%the following appears to be faster than just multiplying by the big kronecker product.
update = zeros(m*r,1);
    for j = 1:r
        update((j-1)*m+1:j*m) =  L*z((j-1)*m+1:j*m);
    end
end
%% check inputs and set defaults.

function [opts,X,A_w,D_w,L_w,r,sz] = processInputs(opts, X, A_w, r)
        [m,n] = size(X);
        sz = [m,n];
        %% Check opts and set defaults.
        % checking valid options are passed.
        valid_options = {'fasta','maxIters','tol','verbose',...
        'earlyStop',...
        'initialization',...
        'randInit','smoothInit',...
        'completion','L',...
        'normalizedLaplacian','LQF',...
        'l2','l1','netnmfscStop'};
        fn = fieldnames(opts);
        for i = 1:length(fn)
            assert(ismember(fn{i},valid_options),...
            ['Invalid option for boxR2R (',fn{i},').',...
            '  Valid choices are: ',...
              'fasta','maxIters',...
              'tol','verbose',...
              'earlyStop',...
              'initialization',...
              'randInit','smoothInit', ...
              'completion','L',  ...
              'normalizedLaplacian', 'LQF', ...
              'l2','l1','netnmfscStop']);
        end
        
        % fasta: inputs to fasta solver.
        if ~isfield(opts,'fasta')
            opts.fasta = struct();
        end
        %  maxIters: The maximum number of outer iterations
        if ~isfield(opts,'maxIters')
            opts.maxIters = 20;
        end
        if ~isfield(opts,'netnmfscStop')
            opts.netnmfscStop = 1e-2;
        end
        % earlyStop: check relative decreases in residuals and gradient
        % norm.
        if ~isfield(opts,'earlyStop')
            opts.earlyStop = true;
        end
        % tol:  The relative decrease in the residuals before the method stops
        if ~isfield(opts,'tol') % Stopping tolerance
            opts.tol = 1e-3;
        end
        % verbose: 'true' => print status information on every iteration
        if ~isfield(opts,'verbose')   
            opts.verbose = false;
        end
        if ~isfield(opts.fasta,'verbose')
            opts.fasta.verbose = opts.verbose;
        end

        

        % initialization.W & .H: matrices to initialize from
        if ~isfield(opts,'initialization')
            opts.initialization = false;
        end
        % randInit: 'true'=> initialize from abs(i.i.d. normal)
        % 'false' => init from nnmf(X,r)

        if ~isfield(opts,'randInit')
            opts.randInit = false;
        end
        % completion: 'true'=>solve problem penalized only by nonzero entries.
        if ~isfield(opts,'completion')   
            opts.completion = true;
        end
        % L: 'true'=> A_w = L
        if ~isfield(opts,'L')   
            opts.L = false;
        end

        % normalizedL: 'true'=> use a normalized Laplacian
        if ~isfield(opts,'normalizedLaplacian')   
            opts.normalizedLaplacian = false;
        end
        %check validity of opts.normalizedLaplacian.
        assert(any([any(opts.normalizedLaplacian==[0,1])...
            islogical(opts.normalizedLaplacian)]),...
            'normalizedLaplacian must be [0,1] or logical.')
        % smoothInit: 'true'=> apply smoothing to first W.
        if ~isfield(opts, 'smoothInit')
            opts.smoothInit = false;
        end
        % LQF: 'isnumeric(opts.LQF)'=> opts.LQF = lambda >=0, graph smoothness penalty
        if ~isfield(opts,'LQF')
            opts.LQF = 0;
        else
            assert(all([isnumeric(opts.LQF) size(opts.LQF)==1 opts.LQF>=0]), ...
            'opts.LQF must be a non-negative number')
        end
        
        % l2: 'isnumeric(opts.l2)'=> opts.l2 = theta >=0, l2 penalty
        if ~isfield(opts,'l2')
            opts.l2 = 0;
        else
            assert(all([isnumeric(opts.l2) size(opts.l2)==1 opts.l2>=0]), ...
            'opts.l2 must be a non-negative number')
        end
        
        % l1: 'isnumeric(opts.l1)'=> opts.l1 = mu >=0, l1 penalty
        if ~isfield(opts,'l1')
            opts.l1 = 0;
%         else
%             assert([isnumeric(opts.l1) size(opts.l1)==1 opts.l1>=0], ...
%             'opts.l1 must be a non-negative number')
        end
        
        % If all regs are zero, make ell2 nonzero in order to choose least
        % norm solutions for stability
        if all([opts.LQF,opts.l2,opts.l1]==0)
            if opts.verbose
                disp('\n Setting l2 penalty to 1e-3 to stabilize unregularized problem')
            end
                opts.l2 = 1e-3;
        end
        
        
        %% check inputs
        if ~issparse(X)
            X = sparse(X);
        end
        if isempty(r)
            r = min([m,n]);
        end
        if opts.L 
            % in this case, L_w was passed by A_w
            L_w = A_w;
            % we set A_w to be the 
            % corresponding adjacency matrix.
            D_w = diag(diag(L_w));
            A_w = D_w-L_w;
        else
            %in this case, L did not get set from opts.L==1
            D_w = sum(A_w,1);
            if ~opts.normalizedLaplacian
                L_w = diag(D_w) - A_w;
            else
                error('normalized Laplacian not currently supported')
            end
        end
  
        if ~issparse(A_w)
            A_w = sparse(A_w);
        end
        if ~issparse(L_w)
            L_w = sparse(L_w);
        end
        assert(r<=min([m,n]), 'input rank r must be leq the smallest dimension of X');
        assert(all(size(A_w) == [m,m]),...
        'Input adjacency A_w must be square and m x m' );
        assert(issparse(A_w),'A_w must be sparse, this error is unreachable');
        assert(issparse(X), 'X must be sparse, this error is unreachable');
end
    