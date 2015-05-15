%Make Fig1 for LOVEFest (Randomer Forest)

close all
clear
clc

fpath = mfilename('fullpath');
findex = strfind(fpath,'/');
rootDir=fpath(1:findex(end-1));
p = genpath(rootDir);
gits=strfind(p,'.git');
colons=strfind(p,':');
for i=0:length(gits)-1
endGit=find(colons>gits(end-i),1);
p(colons(endGit-1):colons(endGit)-1)=[];
end
addpath(p);


Colors = linspecer(6,'sequential');
ScatterColors = {'r' 'b' 'g' 'm'};
Fig_Color = [1 1 1];
LineWidth = 3;
Marker = 'none';
Title = {'(A) Trunk' '(B) Parity' '(C) Multimodal'};
Units = 'pixels';
FigPosition = [0 140 1300 350];
Axis_Left = [75 400 725];
Axis_Bottom = 63;
Axis_Width = 250;
Axis_Height = 250;
Legend_Width = 75;
Legend_Height = 80;
Legend_Left = 1075;
Legend_Bottom = round(Axis_Bottom+Axis_Height/2-Legend_Height/2);
MarkerSize = 24;
Box = 'off';
Box_Legend = 'off';
FontSize = 24;

BasePath = '~/LOVEFest/Figures/fig/';
Filename = {'Trunk_time_vs_d_n100_var1_ntrees1500_ntrials10_v3.fig'...
 'Parity_time_vs_d_n100_ntrees1000_ntrials10_v6.fig'...
 'Multimodal_time_vs_d_n100_var1_ntrees500_ntrials10_v3.fig'};

for i = 1:length(Filename)
    h{i} = openfig(strcat(BasePath,Filename{i}),'invisible');
    grid off
end

h{4} = figure('Visible','On');
set(h{4},'Units','normalized','position',[0 0 1 1]);
set(h{4},'Units','inches');
screenposition = get(h{4},'Position');
set(h{4},...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',screenposition(3:4));

for i = 1:length(h)-1
    ax_old = get(h{i},'CurrentAxes');
    ax_new(i) = subplot(1,3,i);
    axis square
    copyobj(allchild(ax_old),ax_new(i));
    h_lines = allchild(ax_new(i));
    xmax = zeros(1,length(h_lines));
    ymax_ln = zeros(1,length(h_lines));
    for j = 1:length(h_lines)
        set(h_lines(j),'Color',Colors(j,:),'linewidth',LineWidth,'Marker',Marker,'MarkerFaceColor',Colors(j,:),'MarkerEdgeColor',Colors(j,:))
        xmax(j) = max(get(h_lines(j),'XData'));
        ymax_ln(j) = max(get(h_lines(j),'YData'));
    end
    xmax = max(xmax);
    ymax_ax = max(ymax_ln);
    ymax_idx = find(ymax_ln==ymax_ax);
    ymax_ax = ymax_ax + h_lines(ymax_idx).UData(find(h_lines(ymax_idx).YData==ymax_ln(ymax_idx)));
    if i ~= 2
        XTick = logspace(0,log10(xmax),log10(xmax)+1);
        set(ax_new(i),'FontSize',FontSize,'XLim',[0 10^(log10(xmax)+0.1)],'YLim',[0 ymax_ax],'XScale','log','XTick',XTick,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth)
    else
        XTick = 0:2:10;
        set(ax_new(i),'FontSize',FontSize,'XLim',[0 11],'YLim',[0 ymax_ax],'XScale','linear','XTick',XTick,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth)
    end
    title(Title{i})
    xlabel('# Ambient Dim')
    if i == 1
        ylabel('Training Time (sec)')
    end
end
axpos = get(ax_new(end),'Position');
ax_edge = axpos(1) + axpos(3);
hL = legend(ax_new(end),'RF','RerF(d)','RerF(s)','RerF(s+d)','RerF(s+d+r)');
lpos = get(hL,'Position');
set(hL, 'position',[ax_edge+.01 axpos(2)+axpos(4)/2-lpos(4)/2 lpos(3:4)],'Box',Box_Legend,'FontSize',FontSize)

pos_l = get(hL,'Position');
pos_a = get(ax_new(1),'Position');
left_adjust = 1 - (pos_l(1) + pos_l(3));

for i = 1:length(ax_new)
    pos = get(ax_new(i),'Position');
    set(ax_new(i),'Position',[pos(1)+left_adjust pos(2) pos(3) pos(4)])
end

pos = get(hL,'Position');
set(hL,'Position',[pos(1)+left_adjust pos(2) pos(3) pos(4)])

fname = '~/LOVEFest/Figures/Fig2_Time';
save_fig(gcf,fname)
