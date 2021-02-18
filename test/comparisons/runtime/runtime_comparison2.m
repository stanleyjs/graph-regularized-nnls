close all;
clear all;
data_folder = '../../../datasets/';
load([data_folder 'netnmfsc_comparison.mat'])
load('resampled_netnmfsc.mat')

%%
ntrials = 5;
times = zeros(numel(X1),ntrials);
r=5;
for ix = 1:5
    disp(['current multiple ', num2str(ix)])
    clear opts %reset to defaults
    
    % set stable options
    opts.LQF = 10;
    opts.l1 = 0;
    opts.randInit=true;
    opts.fasta.recordObjective=false;
    opts.smoothInit=true;
    opts.fasta.tol=1e-3;
    opts.fasta.maxIters=10;
    opts.maxIters=50;
    opts.fasta.verbose=false;
    %construct initializations
    opts.fasta.adaptive=true;
    opts.fasta.accelerate=true;
    for trial = 1:ntrials
        
        tt = tic;
        z = library_size_normalization(X1{ix});
        [Y,W,H,obj] = boxR2RNNGLS(z,A,r,opts);
        times(ix,trial) = toc(tt);
        disp(obj.obj)
        disp((['%%%% trial completed in ' num2str(times(ix,trial)) ' seconds']))
    end

end
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



