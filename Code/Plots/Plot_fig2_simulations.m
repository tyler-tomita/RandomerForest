close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

LineWidth = 2;
MarkerSize = 12;
FontSize = .2;
TitleFontSize = 18;
axWidth = 1.75;
axHeight = 1.75;
axLeft = [FontSize*4 FontSize*8+axWidth FontSize*12+axWidth*2 FontSize*4,...
    FontSize*8+axWidth FontSize*12+axWidth*2];
axBottom = [FontSize*9+axHeight FontSize*9+axHeight FontSize*9+axHeight,...
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
l.FontSize = 12;
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(1) axBottom(1) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [-2 3];
ax.YLim = [-2 3];
ax.XTick = [-2 0 2];
ax.YTick = [-2 0 2];
ax.XTickLabel = {'-2';'0';'2'};
ax.YTickLabel = {'-2';'0';'2'};

ax = subplot(2,3,2);

for i = 1:length(classifiers)
    cl = classifiers{i};
    if ~strcmp(cl,'rerfdn')
        [minLhat.(cl),minLhatIdx.(cl)] = min(meanLhat.(cl)(end,:,:),[],2);
        for j = 1:length(dims)
            minSemLhat.(cl)(j) = semLhat.(cl)(end,minLhatIdx.(cl)(j),j);
        end
        if i ~= 1
            hLhat(i) = errorbar(dims,minLhat.(cl)(:)'-minLhat.rf(:)',sqrt(minSemLhat.(cl).^2+minSemLhat.rf.^2),'LineWidth',LineWidth);
        end
        hold on
    end
end

if runSims
    Sparse_parity_bayes_error
else
    load Sparse_parity_bayes_error.mat
end

% errorbar(dims,bayes_error,sem_bayes_error,'LineWidth',LineWidth)

title('(B) Sparse Parity')
xlabel('d')
ylabel('Relative Error')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(2) axBottom(2) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [0 105];
% ax.YLim = [0 .5];
l = legend('RerF','RotRF');
l.Location = 'northwest';
l.Box = 'off';

ax = subplot(2,3,3);

for i = 1:length(classifiers)
    cl = classifiers{i};
    if ~strcmp(cl,'rerfdn')
        trainTime.(cl) = nanmean(meanTrainTime.(cl),2);
        semMeanTrainTime.(cl) = nanmean(semTrainTime.(cl),2);
        hTrainTime(i) = errorbar(dims,trainTime.(cl)(:)',semMeanTrainTime.(cl)(:)','LineWidth',LineWidth);
        hold on
    end
end

title('(C) Sparse Parity')
xlabel('d')
ylabel('Train Time (s)')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(3) axBottom(3) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [0 105];
l = legend('RF','RerF','RotRF');
l.Location = 'northwest';
l.Box = 'off';
l.FontSize = 12;

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
ax.XLim = [-5 5];
ax.YLim = [-5 5];
ax.XTick = [-5 0 5];
%ax.YTick = [-5 0 5];
ax.XTickLabel = {'-5';'0';'5'};
%ax.YTickLabel = {'-5';'0';'5'};

ax = subplot(2,3,5);

for i = 1:length(classifiers)
    cl = classifiers{i};
    if ~strcmp(cl,'rerfdn')
        [minLhat.(cl),minLhatIdx.(cl)] = min(meanLhat.(cl)(end,:,:),[],2);
        for j = 1:length(dims)
            minSemLhat.(cl)(j) = semLhat.(cl)(end,minLhatIdx.(cl)(j),j);
        end
        if i ~= 1
            hLhat(i) = errorbar(dims,minLhat.(cl)(:)'-minLhat.rf(:)',sqrt(minSemLhat.(cl).^2+minSemLhat.rf.^2),'LineWidth',LineWidth);
        end
        hold on
    end
end

if runSims
    Trunk_bayes_error
else
    load Trunk_bayes_error.mat
end

% plot(dims,bayes_error,'LineWidth',LineWidth)

title('(E) Trunk')
xlabel('d')
ylabel('Relative Error')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(5) axBottom(5) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [1 550];
% ax.YLim = [0 .15];
ax.XScale = 'log';
ax.XTick = [logspace(0,2,3) 500];
ax.XTickLabel = {'1';'10';'100';'500'};

ax = subplot(2,3,6);

for i = 1:length(classifiers)
    cl = classifiers{i};
    if ~strcmp(cl,'rerfdn')
        trainTime.(cl) = nanmean(meanTrainTime.(cl),2);
        semMeanTrainTime.(cl) = nanmean(semTrainTime.(cl),2);
        hTrainTime(i) = errorbar(dims,trainTime.(cl)(:)',semMeanTrainTime.(cl)(:)','LineWidth',LineWidth);
        hold on
    end
end

%Plot dummy line for bayes just so that it's in the legend
% plot([0 0 0],[0 0 0],'LineWidth',LineWidth,'Visible','off')

title('(F) Trunk')
xlabel('d')
ylabel('Train Time (s)')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(6) axBottom(6) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [1 550];
ax.XScale = 'log';
ax.XTick = [logspace(0,2,3) 500];
ax.XTickLabel = {'1';'10';'100';'500'};

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Fig2_simulations'])