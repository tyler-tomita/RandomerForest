%% Plot Sparse Parity and Trunk

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')
Colors.rf = ColorMap(1,:);
Colors.rerf = ColorMap(3,:);
Colors.frc= ColorMap(9,:);
Colors.rr_rf = ColorMap(11,:);
Colors.xgb = 'k';
LineStyles.rf = '-';
LineStyles.rerf = '-';
LineStyles.frc = '-';
LineStyles.rr_rf = '-';
LineStyles.xgb = '-';
LineWidth = 2;

%% Plot Sparse Parity

load Sparse_parity

for i = 1:length(TestError)
    Classifiers = {'rf','rerf','frc','rr_rf','xgb'};
    for j = 1:length(Classifiers)
        if strcmp(Classifiers{j},'xgb')
            PlotError.(Classifiers{j}) = dlmread('/Users/tyler/Documents/R/Results/Sparse_parity_Lhat.dat');
        else
            ntrials = length(TestError{i}.(Classifiers{j}).Untransformed);
            for trial = 1:ntrials
                PlotError.(Classifiers{j})(trial,i) = ...
                    TestError{i}.(Classifiers{j}).Untransformed(trial);
            end
        end
    end
end

for i = 1:length(Classifiers)
    cl = Classifiers{i};
    hTestError(i) = errorbar(dims,mean(PlotError.(cl)),...
        std(PlotError.(cl))/sqrt(size(PlotError.(cl),1)),...
        'LineWidth',LineWidth,'Color',Colors.(cl),...
        'LineStyle',LineStyles.(cl));
    hold on
end

ax = gca;
ax.LineWidth = LineWidth;
ax.FontSize = 14;
ax.Box = 'off';
xlabel('p')
ylabel('Error Rate')
title('Sparse Parity')

h_leg = legend('RF','RerF','F-RC','RR-RF','XGBoost');
h_leg.Box = 'off';
h_leg.Location = 'southeast';
h_leg.FontSize = 14;

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Fig2_sparse_parity_with_xgboost'])