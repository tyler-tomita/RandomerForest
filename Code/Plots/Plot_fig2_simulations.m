close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

LineWidth = 2;
MarkerSize = 12;
FontSize = .2;
axWidth = 1.75;
axHeight = 1.75;
axLeft = [FontSize*4 FontSize*8+axWidth FontSize*12+axWidth*2 FontSize*4,...
    FontSize*8+axWidth FontSize*12+axWidth*2];
axBottom = [FontSize*8+axHeight FontSize*8+axHeight FontSize*8+axHeight,...
    FontSize*4 FontSize*4 FontSize*4];
figWidth = axLeft(end) + axWidth + FontSize*4;
figHeight = axBottom(1) + axHeight + FontSize*4;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];

runSims = false;

if runSims
    run_Sparse_parity
else
    load Sparse_parity.mat
end

classifiers = fieldnames(meanLhat);

ax = subplot(2,3,1);

n = 1000;
d = 2;
dgood = min(3,d);
Sigma = 1/32*ones(1,d);
Mu = sparse(n,d);

for jj = 1:n
    Mu(jj,:) = binornd(1,0.5,1,d);
    X(jj,1:d) = mvnrnd(Mu(jj,:),Sigma);
end

nones = sum(Mu(:,1:dgood),2);
Y = mod(nones,2);

plot(X(Y==0,1),X(Y==0,2),'b.',X(Y==1,1),X(Y==1,2),'r.','MarkerSize',MarkerSize)

title('(A) Sparse Parity')
xlabel('X_1')
ylabel('X_2')
l = legend('Class 1','Class 2');
l.Location = 'northwest';
l.Box = 'off';
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(1) axBottom(1) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [-2 3];
ax.YLim = [-2 3];

ax = subplot(2,3,2);

for i = 1:length(classifiers)
    cl = classifiers{i};
    [minLhat.(cl),minLhatIdx.(cl)] = min(meanLhat.(cl)(end,:,:),[],2);
    for j = 1:length(dims)
        minSemLhat.(cl)(j) = semLhat.(cl)(end,minLhatIdx.(cl)(j),j);
    end
    hLhat(i) = errorbar(dims,minLhat.(cl)(:)',minSemLhat.(cl),'LineWidth',LineWidth);
    hold on
end

title('(B) Sparse Parity')
xlabel('d')
ylabel('Error Rate')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(2) axBottom(2) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [0 105];
ax.YLim = [0 .5];

ax = subplot(2,3,3);

for i = 1:length(classifiers)
    cl = classifiers{i};
    trainTime.(cl) = nanmean(meanTrainTime.(cl),2);
    semMeanTrainTime.(cl) = nanmean(semTrainTime.(cl),2);
    hTrainTime(i) = errorbar(dims,trainTime.(cl)(:)',semMeanTrainTime.(cl)(:)','LineWidth',LineWidth);
    hold on
end

title('(C) Sparse Parity')
xlabel('d')
ylabel('Training Time (sec)')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(3) axBottom(3) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [0 105];

clear hLhat hTrainTime minLhat minSemLhat trainTime semMeanTrainTime

if runSims
    run_Trunk
else
    load Trunk.mat
end

classifiers = fieldnames(meanLhat);

ax = subplot(2,3,4);

n = 100;
d = 2;

d_idx = 1:d;
mu1 = 1./sqrt(d_idx);
mu0 = -1*mu1;
Mu = cat(1,mu0,mu1);
Sigma = ones(1,d);
obj = gmdistribution(Mu,Sigma);
[X,idx] = random(obj,n);
Class = [0;1];
Y = Class(idx);

plot(X(Y==0,1),X(Y==0,2),'b.',X(Y==1,1),X(Y==1,2),'r.','MarkerSize',MarkerSize)

title('(D) Trunk')
xlabel('X_1')
ylabel('X_2')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(4) axBottom(4) axWidth axHeight];
ax.Box = 'off';

ax = subplot(2,3,5);

for i = 1:length(classifiers)
    cl = classifiers{i};
    [minLhat.(cl),minLhatIdx.(cl)] = min(meanLhat.(cl)(end,:,:),[],2);
    for j = 1:length(dims)
        minSemLhat.(cl)(j) = semLhat.(cl)(end,minLhatIdx.(cl)(j),j);
    end
    hLhat(i) = errorbar(dims,minLhat.(cl)(:)',minSemLhat.(cl),'LineWidth',LineWidth);
    hold on
end

title('(E) Trunk')
xlabel('d')
ylabel('Error Rate')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(5) axBottom(5) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [1 550];
ax.YLim = [0 .15];
ax.XScale = 'log';
ax.XTick = [logspace(0,2,3) 500];

ax = subplot(2,3,6);

for i = 1:length(classifiers)
    cl = classifiers{i};
    trainTime.(cl) = nanmean(meanTrainTime.(cl),2);
    semMeanTrainTime.(cl) = nanmean(semTrainTime.(cl),2);
    hTrainTime(i) = errorbar(dims,trainTime.(cl)(:)',semMeanTrainTime.(cl)(:)','LineWidth',LineWidth);
    hold on
end

title('(F) Trunk')
xlabel('d')
ylabel('Training Time (sec)')
l = legend('RF','RerF','RerFd','Rotation RF');
l.Location = 'northwest';
l.Box = 'off';
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(6) axBottom(6) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [1 550];
ax.XScale = 'log';
ax.XTick = [logspace(0,2,3) 500];

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Fig2_simulations'])