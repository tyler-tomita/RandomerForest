close all
clear
clc

open ROFLMAO_fig2_posteriors.fig

load('purple2green')
ColorMap2 = interpolate_colormap(ColorMap,64,true);
ColorMap3 = jet;
ColorMap3 = [ColorMap3(1,:);ones(1,3);ColorMap3(end,:)];
ColorMap3 = interpolate_colormap(ColorMap3,64,true);

Colors.rf = ColorMap(1,:);
Colors.rfr = ColorMap(1,:);
Colors.frc= ColorMap(9,:);
Colors.frcr = ColorMap(9,:);
Colors.rr_rf = ColorMap(3,:);
Colors.rr_rfr = ColorMap(3,:);
LineStyles.rf = '-';
LineStyles.rfr = ':';
LineStyles.frc = '-';
LineStyles.frcr = ':';
LineStyles.rr_rf = '-';
LineStyles.rr_rfr = ':';

LineWidth = 2;
FontSize = .15;
axWidth = .75;
axHeight = .75;
cbWidth = .1;
cbHeight = axHeight;
axLeft = repmat([FontSize*3,FontSize*5+axWidth,FontSize*7+2*axWidth],1,8);
% axLeft = [FontSize*3,FontSize*5+axWidth,FontSize*7+2*axWidth,FontSize*3,...
%     FontSize*5+axWidth,FontSize*7+2*axWidth,FontSize*3,...
%     FontSize*5+axWidth,FontSize*7+2*axWidth,FontSize*3,...
%     FontSize*5+axWidth,FontSize*7+2*axWidth,FontSize*3,...
%     FontSize*5+axWidth,FontSize*7+2*axWidth];
axBottom = [(FontSize*13.5+axHeight*8)*ones(1,3),(FontSize*12+axHeight*7)*ones(1,3),(FontSize*10.5+axHeight*6)*ones(1,3),(FontSize*9+axHeight*5)*ones(1,3),...
    (FontSize*7.5+axHeight*4)*ones(1,3),(FontSize*6+axHeight*3)*ones(1,3),...
    (FontSize*4.5+axHeight*2)*ones(1,3),(FontSize*3+axHeight)*ones(1,3)];
axHeight2 = axHeight;
axWidth2 = axLeft(end) - axLeft(1);
axLeft2 = axLeft(1) + axWidth(1)/2;
axBottom2 = FontSize*1.5;
cbLeft = axLeft + axWidth + FontSize/2;
cbBottom = axBottom;

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

ax = findobj(gcf,'type','axes');
for i = 2:16
    colormap(ax(i),ColorMap2)
end

for i = 17:19
    colormap(ax(i),ColorMap3)
end

ln = findobj(ax(1),'type','errorbar');

for i = 1:length(ln)
    if strcmp(ln(i).DisplayName,'RF')
        ln(i).Color = Colors.rf;
        ln(i).LineStyle = LineStyles.rf;
    elseif strcmp(ln(i).DisplayName,'F-RC')
        ln(i).Color = Colors.frc;
        ln(i).LineStyle = LineStyles.frc;
    elseif strcmp(ln(i).DisplayName,'Frank')
        ln(i).Color = Colors.frcr;
        ln(i).LineStyle = LineStyles.frcr;
    elseif strcmp(ln(i).DisplayName,'RR-RF')
        ln(i).Color = Colors.rr_rf;
        ln(i).LineStyle = LineStyles.rr_rf;
    elseif strcmp(ln(i).DisplayName,'RR-RF(r)')
        ln(i).Color = Colors.rr_rfr;
        ln(i).LineStyle = LineStyles.rr_rfr;
    end
end

legWidth = ax(1).Position(3)/2;
legHeight = ax(1).Position(4);
legLeft = ax(1).Position(1) + ax(1).Position(3) - FontSize*2;
legBottom = ax(1).Position(2);

[lh,objh] = legend(ax(1),'RF','F-RC','Frank','RR-RF','RR-RF(r)');
lh.Units = 'inches';
lh.Position = [legLeft legBottom legWidth legHeight];
lh.Box = 'off';

for i = 6:length(objh)
    objh(i).Children.Children(2).XData = [(objh(i).Children.Children(2).XData(2)-objh(i).Children.Children(2).XData(1))*0.5+objh(i).Children.Children(2).XData(1),objh(i).Children.Children(2).XData(2)];
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/ROFLMAO_fig2_posteriors_2017_01_23'],{'fig','pdf','png'})