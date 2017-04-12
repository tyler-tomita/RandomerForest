close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')
Colors.rerf= ColorMap(10,:);
Colors.rr_rf = ColorMap(4,:);

LineWidth = 2;
MarkerSize = 12;
FontSize = .2;
axWidth = 2;
axHeight = 1.5;
axLeft = repmat([FontSize*4,FontSize*7+axWidth],1,4);
axBottom = [(FontSize*12+axHeight*3)*ones(1,2)...
    (FontSize*9+axHeight*2)*ones(1,2),(FontSize*6+axHeight)*ones(1,2),...
    FontSize*3*ones(1,2)];
legWidth = axWidth;
legHeight = axHeight;
legLeft = axLeft(end) + axWidth - FontSize;
legBottom = mean([axBottom(4) axBottom(6)]);
figWidth = legLeft + legWidth;
figHeight = axBottom(1) + axHeight + FontSize*2;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

ps = 2.^(1:6);
ntrials = 1e3;
deltas = [1,10,22.5,45] + 0.1; %add 0.1 b/c of slight imprecision in angle calculation method

k = 1;

isClose1 = false(length(deltas),length(ps),ntrials);
isClose2 = false(length(deltas),length(ps),ntrials);
isClose3 = false(length(deltas),length(ps),ntrials);

for i = 1:length(deltas)
    delta = deltas(i);
    for j = 1:length(ps)
        p = ps(j);
        target_vector = [1 zeros(1,p-1)];
        for trial = 1:ntrials
            R1 = random_rotation(p);
            R2 = randmat(p,p,'binary',1/p);
            R3 = randmat(p,p^2,'binary',1/p);
            ncol = size(R2,2);
            ncol2= size(R3,2);
            success1 = false(p,1);
            success2 = false(p,1);
            success3 = false(p,1);
            for c = 1:p
                theta1 = vector_angle(target_vector,R1(:,c),'degrees');
                success1(c) = (theta1 <= delta) | (theta1 >= 180 - delta);
                if c <= ncol
                    theta2 = vector_angle(target_vector,R2(:,c),'degrees');
                    success2(c) = (theta2 <= delta) | (theta2 >= 180 - delta);
                end
                if c <= ncol2
                    theta3 = vector_angle(target_vector,R3(:,c),'degrees');
                    success3(c) = (theta3 <= delta) | (theta3 >= 180 - delta);
                end
            end
            isClose1(i,j,trial) = any(success1);
            isClose2(i,j,trial) = any(success2);
            isClose3(i,j,trial) = any(success3);
        end
    end
    
    ax((i-1)*2+k) = axes;
    hold on
    errorbar(ps,mean(isClose2(i,:,:),3),std(isClose2(i,:,:),0,3)/sqrt(ntrials),...
        'LineStyle','-','Color',Colors.rerf,'LineWidth',2)
    errorbar(ps,mean(isClose3(i,:,:),3),std(isClose3(i,:,:),0,3)/sqrt(ntrials),...
        'LineStyle',':','Color',Colors.rerf,'LineWidth',2)
    errorbar(ps,mean(isClose1(i,:,:),3),std(isClose1(i,:,:),0,3)/sqrt(ntrials),...
        'LineStyle','-','Color',Colors.rr_rf,'LineWidth',2)
    
    xlabel('p')
    ylabel(['P(\Theta \leq ' num2str(delta-0.1) '°)'])
    if i == 1
        title('v^* = (1,0,...0)')
    end

    ax((i-1)*2+k).LineWidth = LineWidth;
    ax((i-1)*2+k).FontUnits = 'inches';
    ax((i-1)*2+k).FontSize = FontSize;
    ax((i-1)*2+k).Units = 'inches';
    ax((i-1)*2+k).Position = [axLeft((i-1)*2+k) axBottom((i-1)*2+k) axWidth axHeight];
    ax((i-1)*2+k).Box = 'off';
    YT = ax((i-1)*2+k).YTick;
    ax((i-1)*2+k).YLim = [ax((i-1)*2+k).YLim(1)-diff(ax((i-1)*2+k).YLim)*0.1 ax((i-1)*2+k).YLim(2)];
    ax((i-1)*2+k).YTick = YT;
    ax((i-1)*2+k).XScale = 'log';
%     ax((i-1)*2+k).XLim = [10^(log10(min(ns{j}))-0.1) 10^(log10(max(ns{j}))+0.1)];
%     ax((i-1)*2+k).XScale = 'log';
%     ax((i-1)*2+k).XTick = ns{j};
%     ax((i-1)*2+k).XTickLabel = cellstr(num2str(ns{j}'))';
%     ax((i-1)*2+k).XTickLabelRotation = 45;
%     ax((i-1)*2+k).YLim = [min(ErrorMatrix(:)) max(ErrorMatrix(:))];
end

k = 2;

isClose1 = false(length(deltas),length(ps),ntrials);
isClose2 = false(length(deltas),length(ps),ntrials);
isClose3 = false(length(deltas),length(ps),ntrials);

for i = 1:length(deltas)
    delta = deltas(i);
    for j = 1:length(ps)
        p = ps(j);
        target_vector = ones(1,p);
        for trial = 1:ntrials
            R1 = random_rotation(p);
            R2 = randmat(p,p,'binary',1/p);
            R3 = randmat(p,p^2,'binary',1/p);
            ncol = size(R2,2);
            ncol2= size(R3,2);
            success1 = false(p,1);
            success2 = false(p,1);
            success3 = false(p,1);
            for c = 1:p
                theta1 = vector_angle(target_vector,R1(:,c),'degrees');
                success1(c) = (theta1 <= delta) | (theta1 >= 180 - delta);
                if c <= ncol
                    theta2 = vector_angle(target_vector,R2(:,c),'degrees');
                    success2(c) = (theta2 <= delta) | (theta2 >= 180 - delta);
                end
                if c <= ncol2
                    theta3 = vector_angle(target_vector,R3(:,c),'degrees');
                    success3(c) = (theta3 <= delta) | (theta3 >= 180 - delta);
                end
            end
            isClose1(i,j,trial) = any(success1);
            isClose2(i,j,trial) = any(success2);
            isClose3(i,j,trial) = any(success3);
        end
    end
    
    ax((i-1)*2+k) = axes;
    hold on
    errorbar(ps,mean(isClose2(i,:,:),3),std(isClose2(i,:,:),0,3)/sqrt(ntrials),...
        'LineStyle','-','Color',Colors.rerf,'LineWidth',2)
    errorbar(ps,mean(isClose3(i,:,:),3),std(isClose3(i,:,:),0,3)/sqrt(ntrials),...
        'LineStyle',':','Color',Colors.rerf,'LineWidth',2)
    errorbar(ps,mean(isClose1(i,:,:),3),std(isClose1(i,:,:),0,3)/sqrt(ntrials),...
        'LineStyle','-','Color',Colors.rr_rf,'LineWidth',2)
    
    xlabel('p')
%     ylabel(['P(\Theta \leq ' num2str(delta-0.1) '°)'])
    if i == 1
        title('v^* = (1,1,...1)')
    end

    ax((i-1)*2+k).LineWidth = LineWidth;
    ax((i-1)*2+k).FontUnits = 'inches';
    ax((i-1)*2+k).FontSize = FontSize;
    ax((i-1)*2+k).Units = 'inches';
    ax((i-1)*2+k).Position = [axLeft((i-1)*2+k) axBottom((i-1)*2+k) axWidth axHeight];
    ax((i-1)*2+k).Box = 'off';
    YT = ax((i-1)*2+k).YTick;
    ax((i-1)*2+k).YLim = [ax((i-1)*2+k).YLim(1)-diff(ax((i-1)*2+k).YLim)*0.1 ax((i-1)*2+k).YLim(2)];
    ax((i-1)*2+k).YTick = YT;
    ax((i-1)*2+k).XScale = 'log';
%     ax((i-1)*2+k).XLim = [10^(log10(min(ns{j}))-0.1) 10^(log10(max(ns{j}))+0.1)];
%     ax((i-1)*2+k).XScale = 'log';
%     ax((i-1)*2+k).XTick = ns{j};
%     ax((i-1)*2+k).XTickLabel = cellstr(num2str(ns{j}'))';
%     ax((i-1)*2+k).XTickLabelRotation = 45;
%     ax((i-1)*2+k).YLim = [min(ErrorMatrix(:)) max(ErrorMatrix(:))];
end

[lh,objh] = legend('RerF (d = p)','RerF (d = p^2)','RR-RF');
lh.Box = 'off';
lh.Units = 'inches';
lh.Position = [legLeft legBottom legWidth legHeight];
for j = length(objh)/2+1:length(objh)
    objh(j).Children.Children(2).XData = [(objh(j).Children.Children(2).XData(2)-objh(j).Children.Children(2).XData(1))*.75+objh(j).Children.Children(2).XData(1),objh(j).Children.Children(2).XData(2)];
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/PAMI_fig1_split_probability'],{'fig','pdf','png'})