close all
clear
clc

load('purple2green')
LineWidth = 2;
FontSize = 14;

d = dir('*.mat');
L = zeros(length(d),5);
fpr = zeros(length(d),5);
ntest = zeros(length(d),1);

for i = 1:length(d) 
    Fold = regexp(d(i).name,'\d*','match');
    Fold = Fold{1};
    
    load(d(i).name)

    L(i,1) = TestError.rf;
    L(i,2) = TestError.rerf;
    L(i,3) = TestError.rr_rf;

    fpr(i,1) = ConfusionMatrix.rf(1,2)/sum(ConfusionMatrix.rf(:));
    fpr(i,2) = ConfusionMatrix.rerf(1,2)/sum(ConfusionMatrix.rerf(:));
    fpr(i,3) = ConfusionMatrix.rr_rf(1,2)/sum(ConfusionMatrix.rr_rf(:));
    
    load(['../2017.01.29/Fold_' Fold '_dense.mat'])
    
    L(i,4) = TestError.rerf_dense;
    
    fpr(i,4) = ConfusionMatrix.rerf_dense(1,2)/sum(ConfusionMatrix.rerf_dense(:));
    
    load(['../2017.01.30/Fold_' Fold '_dense.mat'])
    
    L(i,5) = TestError.rerfc;
    
    fpr(i,5) = ConfusionMatrix.rerfc(1,2)/sum(ConfusionMatrix.rerfc(:));
    
    X = dlmread(['../../Data/extracted/' Fold '.test.dat']);
    ntest(i) = size(X,1);
end

L = L(:,[1,2,4,5,3]);
fpr = fpr(:,[1,2,4,5,3]);

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

%%%%%%%

% d = dir('*.mat');
% L = zeros(length(d),1);
% fpr = zeros(length(d),1);
% ntest = zeros(length(d),1);
% 
% for i = 1:length(d) 
%     S1 = load(d(i).name);
%     Fold = regexp(d(i).name,'\d*','match');
%     Fold = Fold{1};
%     S2 = load(['../2017.01.27/Fold_' Fold '.mat']);
%     S = S2;
%     S.TestError.rerf_dense = S1.TestError.rerf_dense;
%     
%     L(i) = TestError.rerf_dense;
% 
%     fpr(i) = ConfusionMatrix.rerf_dense(1,2)/sum(ConfusionMatrix.rerf_dense(:));
%     
%     X = dlmread(['../../Data/extracted/' Fold '.test.dat']);
%     ntest(i) = size(X,1);
% end
% 
% Weights = ntest/sum(ntest);
% 
% MeanTestError = sum(L.*Weights);
% Meanfpr = sum(fpr.*Weights);