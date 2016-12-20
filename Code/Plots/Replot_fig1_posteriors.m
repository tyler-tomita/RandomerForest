%% Plot posterior heat maps

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('red2blue')
cmap.Posteriors = interpolate_colormap(ColorMap(2:end-1,:),64,true);

load('purple2green')
cmap.Phats = interpolate_colormap(ColorMap(2:end-1,:),64,true);

cmap.Algorithms = ColorMap([1,3,9,11],:);

open Fig1_posteriors.fig

ax = flipud(findobj(gcf,'Type','Axes'));

for i = 1:3
    colormap(ax(i),cmap.Posteriors);
end

for i = 4:length(ax)-1
    colormap(ax(i),cmap.Phats);
end

hlines = findobj(ax(end),'Type','ErrorBar');

for i = 1:length(hlines)
    if strcmp(hlines(i).DisplayName,'RF')
        hlines(i).Color = cmap.Algorithms(1,:);
    elseif strcmp(hlines(i).DisplayName,'RerF') || strcmp(hlines(i).DisplayName,'RerF(r)')
        hlines(i).Color = cmap.Algorithms(2,:);
    elseif strcmp(hlines(i).DisplayName,'RerF2') || strcmp(hlines(i).DisplayName,'RerF2(r)')
        hlines(i).Color = cmap.Algorithms(3,:);
    else
        hlines(i).Color = cmap.Algorithms(4,:);
    end
end

hleg_old = findobj(gcf,'Type','Legend');
u = hleg_old.Units;
po = hleg_old.Position;
fs = hleg_old.FontSize;
[hleg_new,objh] = legend(ax(end),'RF','RerF','RerF(r)','RerF2','RerF2(r)','F-RC','Frank');
hleg_new.Units = u;
hleg_new.Position = po;
hleg_new.FontSize = fs;
hleg_new.Box = 'off';

for i = 8:length(objh)
    objh(i).Children.Children(2).XData = [(objh(i).Children.Children(2).XData(2)-objh(i).Children.Children(2).XData(1))*.75+objh(i).Children.Children(2).XData(1),objh(i).Children.Children(2).XData(2)];
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Fig1_posteriors_v2'])