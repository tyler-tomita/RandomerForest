clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

C = [0 1 1;0 1 0;1 0 1;1 0 0;0 0 0;1 .5 0];
Colors.rf = C(1,:);
Colors.rerf = C(2,:);
Colors.rf_rot = C(3,:);
Colors.rerfr = C(4,:);
Colors.frc = C(5,:);
Colors.rerfdn = C(6,:);
MarkerSize = 12;
FontSize = .2;
TitleFontSize = 18;
axWidth = 1.75;
axHeight = 1.75;
axLeft = FontSize*4;
axBottom = FontSize*4;
figWidth = axLeft(end) + axWidth + FontSize*4;
figHeight = axBottom(1) + axHeight + FontSize*4;

% fig = figure;
% fig.Units = 'inches';
% fig.PaperUnits = 'inches';
% fig.Position = [0 0 figWidth figHeight];
% fig.PaperPosition = [0 0 figWidth figHeight];
% fig.PaperSize = [figWidth figHeight];

Markers = {'none','none','none','none','o'};
LineStyles = {'-','--',':','-.','none'};

ax = axis;
hold on

load Trunk_rerf_frc_parameter_sweep_d50_n1000

mtrys = ceil(d.^[0 1/2 1 1.5 2]);
nmixs = 2:6;

LineWidth = 1:length(mtrys);

for i = 1:length(mtrys)
    err.rerf = Lhat.rerf(3,length(nmixs)*(i-1)+1:length(nmixs)*(i-1)+length(nmixs),:);
    err.frc = Lhat.frc(3,length(nmixs)*(i-1)+1:length(nmixs)*(i-1)+length(nmixs),:);
    errorbar(nmixs,mean(err.rerf,3),std(err.rerf,0,3)/sqrt(size(err.rerf,3)-1),'LineWidth',LineWidth(i),'Color',Colors.rerf)%,'Marker',Markers{i},'LineStyle',LineStyles{i})
    errorbar(nmixs,mean(err.frc,3),std(err.frc,0,3)/sqrt(size(err.frc,3)-1),'LineWidth',LineWidth(i),'Color',Colors.frc)%,'Marker',Markers{i},'LineStyle',LineStyles{i})
    l{2*i-1} = ['RerF (mtry = ' num2str(mtrys(i)) ')'];
    l{2*i} = ['FRC (mtry = ' num2str(mtrys(i)) ')'];
end

xlabel('# linearly combined variables')
ylabel('Error Rate')
title(sprintf('Trunk (n = %d, p = %d)',1000,d))
lg = legend(l);
lg.Box = 'off';
lg.Location = 'eastoutside';
ax = gca;
ax.LineWidth = LineWidth(1);
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
% ax.Units = 'inches';
% %ax.Position = [axLeft axBottom axWidth axHeight];
ax.Box = 'off';

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Trunk_error_vs_sparsity_d50_n1000'])
