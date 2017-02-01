close all
clear
clc

load('purple2green')
LineWidth = 2;
FontSize = 14;

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

resultsPath = [rerfPath 'RandomerForest/Karchin_aneuploidy/Results/2017.01.31/'];
dataPath = [rerfPath 'RandomerForest/Karchin_aneuploidy/Data/extracted/'];
d = dir([resultsPath '*.mat']);

L = zeros(length(d),5);
fpr = zeros(length(d),5);
ntest = zeros(length(d),1);

for i = 1:length(d) 
    Fold = regexp(d(i).name,'\d*','match');
    Fold = Fold{1};
    
    load([resultsPath d(i).name])

    L(i,1) = TestError.rf;
    L(i,2) = TestError.rerf;
    L(i,3) = TestError.rerf_dense;
    L(i,4) = TestError.rerfc;
    L(i,5) = TestError.rr_rf;


    fpr(i,1) = ConfusionMatrix.rf(1,2)/sum(ConfusionMatrix.rf(:));
    fpr(i,2) = ConfusionMatrix.rerf(1,2)/sum(ConfusionMatrix.rerf(:));
    fpr(i,3) = ConfusionMatrix.rerf_dense(1,2)/sum(ConfusionMatrix.rerf_dense(:));
    fpr(i,4) = ConfusionMatrix.rerfc(1,2)/sum(ConfusionMatrix.rerfc(:));
    fpr(i,5) = ConfusionMatrix.rr_rf(1,2)/sum(ConfusionMatrix.rr_rf(:));
    
    Scores = [];
    Scores(:,1) = TestScores.rf(:,2);
    Scores(:,2) = TestScores.rerf(:,2);
    Scores(:,3) = TestScores.rerf_dense(:,2);
    Scores(:,4) = TestScores.rerfc(:,2);
    Scores(:,5) = TestScores.rr_rf(:,2);
    
    filename = fullfile(resultsPath, [Fold '.scores.dat']);
    fid = fopen(filename, 'wt');
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\n','RF','RerF(rho=.0127)','RerF(rho=.0253)','RerF(rho=1)','RR-RF');  % header
    fclose(fid);
    dlmwrite(filename,Scores,'delimiter','\t','precision','%.12f','-append');
    
    X = dlmread([dataPath Fold '.test.dat']);
    ntest(i) = size(X,1);
end

% L = L(:,[1,2,4,5,3]);
% fpr = fpr(:,[1,2,4,5,3]);

Weights = ntest/sum(ntest);
Weights = repmat(Weights,1,size(L,2));

TestError.Mean = sum(L.*Weights);
TestError.SEM = sqrt(var(L,Weights(:,1)))/sqrt(length(d));
FalsePosRate.Mean = sum(fpr.*Weights);
FalsePosRate.SEM = sqrt(var(fpr,Weights(:,1)))/sqrt(length(d));

ax = subplot(2,1,1);
barwitherr(TestError.SEM,TestError.Mean,'LineWidth',LineWidth,...
    'FaceColor',ColorMap(10,:),'EdgeColor',ColorMap(11,:))
title('Aneuploidy Classification')
ax.XTickLabel = {'RF','RerF(rho=.0127)','RerF(rho=.0253)','RerF(rho=1)','RR-RF'};
ax.XTickLabelRotation = 45;
ax.LineWidth = LineWidth;
ax.FontSize = FontSize;
ylabel('21-fold CV Error')
eb = findobj(ax,'Type','ErrorBar');
eb.LineWidth = LineWidth;
eb.Color = ColorMap(2,:);

ax = subplot(2,1,2);
barwitherr(FalsePosRate.SEM,FalsePosRate.Mean,'LineWidth',LineWidth,...
    'FaceColor',ColorMap(10,:),'EdgeColor',ColorMap(11,:))
ax.XTickLabel = {'RF','RerF(rho=.0127)','RerF(rho=.0253)','RerF(rho=1)','RR-RF'};
ax.XTickLabelRotation = 45;
ax.LineWidth = LineWidth;
ylabel('21-fold CV FPR')
eb = findobj(ax,'Type','ErrorBar');
eb.LineWidth = LineWidth;
eb.Color = ColorMap(2,:);
ax.FontSize = FontSize;

save_fig(gcf,'~/RandomerForest/Karchin_aneuploidy/Plots/Error_FPR',{'fig','pdf','png'})