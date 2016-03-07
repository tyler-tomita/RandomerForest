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
LineWidth = 2;
MarkerSize = 12;
FontSize = .2;
TitleFontSize = 18;
axWidth = 1.75;
axHeight = 1.75;
axLeft = FontSize*4;
axBottom = FontSize*4;
figWidth = axLeft(end) + axWidth + FontSize*4;
figHeight = axBottom(1) + axHeight + FontSize*4;

fig = gcf;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

load Trunk_vary_n

mtrys = ceil(d.^[0 1/4 1/2 3/4 1]);

for i = 1:length(mtrys)
    [minLhat.rf,minLhatIdx.rf] = min(mean(Lhat.rf,3),[],2);
    [minLhat.rerf,minLhatIdx.rerf] = min(mean(Lhat.rerf,3),[],2);
    for j = 1:length(ns)
        sem.rf(j) = std(Lhat.rf(j,minLhatIdx.rf(j),:))/sqrt(size(Lhat.rf,3)-1);
        sem.rerf(j) = std(Lhat.rerf(j,minLhatIdx.rerf(j),:))/sqrt(size(Lhat.rerf,3)-1);
    end
    errorbar(ns,minLhat.rf,sem.rf,'LineWidth',LineWidth,'Color',Colors.rf)
    hold on
    errorbar(ns,minLhat.rf,sem.rerf,'LineWidth',LineWidth,'Color',Colors.rerf)
end

load Trunk_bayes_error
bayes_error = bayes_error(end-1);
plot([1 ns(end)],bayes_error*ones(1,2),':k')

ax = gca;

xlabel('n')
ylabel('Error Rate')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
%ax.Position = [axLeft axBottom axWidth axHeight];
ax.Box = 'off';
ax.XLim = [1 ns(end)];
ax.XScale = 'log';

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Trunk_vary_n'])