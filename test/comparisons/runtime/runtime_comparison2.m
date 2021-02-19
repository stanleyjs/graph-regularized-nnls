close all;
clear all;
data_folder = '../../../datasets/';
load([data_folder 'netnmfsc_comparison.mat'])
load('resampled_netnmfsc.mat')

%%
clear opts
ntrials = 3;
times = zeros(numel(X1),ntrials);
r=5;
outputs = {}
for ix = 5
    disp(['current multiple ', num2str(ix)])
    clear opts %reset to defaults
    
    % set stable options
opts.LQF = 10;
opts.l1 = 0;
opts.randInit=false;
opts.fasta.recordObjective=false;
opts.smoothInit=true;
opts.fasta.tol=1e-3;
opts.fasta.maxIters = 100;
opts.maxIters=1000;
opts.earlyStop=false;
opts.fasta.verbose=false;
opts.fasta.adaptive=true;
opts.fasta.accelerate=true;

opts.netnmfscStop=false;
opts.verbose=true;
    for trial = 1:ntrials
        outputs{ix}.obj = {};
        tt = tic;
        z = library_size_normalization(X1{ix});
        [Y,W,H,obj,W0,H0] = boxR2RNNGLS(z,A,r,opts);
        outputs{ix}.W0 = W0;
        outputs{ix}.H0 = H0;
        outputs{ix}.obj{trial} = obj;
        times(ix,trial) = toc(tt);
        disp(obj.obj)
        disp((['%%%% trial completed in ' num2str(times(ix,trial)) ' seconds']))
    end

end


%%
clear opts
opts.LQF = 10;
opts.l1 = 0;
opts.initialization.W = W0;
opts.initialization.H = H0;
opts.randInit=false;
opts.fasta.recordObjective=false;
opts.smoothInit=false;
opts.fasta.tol=1e-3;
opts.fasta.maxIters = 2000;
opts.maxIters=100;
opts.earlyStop=false;
opts.fasta.verbose=true;
opts.fasta.adaptive=true;
opts.fasta.accelerate=false;

opts.fasta.tau=1000;
opts.verbose=true;
tt = tic;
opts.netnmfscStop=false;

[Y,W,H,obj,W0,H0] = boxR2RNNGLS(z,A,r,opts);
t = toc(tt);
%%
save('smooth_initialization_objectives_large','W0','H0','obj');
%%
load('python_runtimes_obj')
%%
figure
elapsed_time = [0 cumsum([ obj.time.iter])];
 objective = obj.obj(2:end);
 elapsed_pytimes = [0 ;(pytimes)];
 plot(elapsed_pytimes(1:100:find(elapsed_pytimes>60,1)),py_obj(1:100:find(elapsed_pytimes>60,1)),'*-','LineWidth',2)
 hold on
plot(elapsed_time(elapsed_time<60),objective(elapsed_time<60),'*-','LineWidth',2)
% 

% hold on
xlabel('time elapsed','FontSize',20)
ylabel('objective function','FontSize',20)

set(gca,'YScale','log')
legend({'netnmfsc','r2r'},'FontSize',20)
title("Non-negative factorization: objective vs. time, 100x1440",'FontSize',20)
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters large_runtime_comparison2

%%
figure
elapsed_time = [0 cumsum([ obj.time.iter])];
 objective = obj.obj(2:end);
 elapsed_pytimes = [0 ;(pytimes)];
 plot(elapsed_pytimes(1:100:end),py_obj(1:100:end),'*-','LineWidth',2)
 hold on
plot(elapsed_time,objective,'*-','LineWidth',2)
% 

% hold on
xlabel('time elapsed','FontSize',20)
ylabel('objective function','FontSize',20)

set(gca,'YScale','log')
legend({'netnmfsc','r2r'},'FontSize',20)
title("Non-negative factorization: objective vs. time, 100x1440",'FontSize',20)
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters large_runtime_comparison2
%%
save('matlab_times', 'times')
%%
close all
clear all;
load matlab_times
load python_times
load resampled_netnmfsc
N = zeros(numel(X1),1);
for ix= 1:numel(X1)
    N(ix) = size(X1{ix},2);
end

matlab_means = mean(times,2);
matlab_std = std(times');
py_means = mean(pytimes,2);
py_std= std(pytimes');
figure
y1=errorbar(N,matlab_means,matlab_std,'r.-','MarkerSize',20,'CapSize',18,'LineWidth',1)
hold on;
plot(N,matlab_means, 'r-', 'LineWidth',2)
for ix = 1:numel(X1)
    scatter(repmat(N(ix),numel(matlab_means),1),times(ix,:),200,'r*')
end



y2=errorbar(N,py_means,py_std,'b.-','MarkerSize',20,'CapSize',18,'LineWidth',1)
hold on;
plot(N,py_means, 'b-', 'LineWidth',2)
for ix = 1:numel(X1)
    scatter(repmat(N(ix),numel(py_means),1),pytimes(ix,:),200,'b*')
end
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

xlabel('Number of columns','FontSize', 24)
ylabel('Time (s)','FontSize', 24)

title('R2R with improvements: still scales poorly','FontSize', 24)
legend([y2,y1],{'netnmfsc','r2r'})
axis('equal')

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters comparison1



