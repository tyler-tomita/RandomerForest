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

Colors = linspecer(5,'sequential');
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
Filename = {'Trunk_time_vs_d_n100_var1_ntrees1000_ntrials10_v2.fig'...
 'Parity_time_vs_d_n100_ntrees1000_ntrials10_v5.fig'...
 'Multimodal_time_vs_d_n100_var1_ntrees10_ntrials10_v2.fig'};

for i = 1:length(Filename)
    h{i} = openfig(strcat(BasePath,Filename{i}),'invisible');
    grid off
end

h{4} = figure('Visible','On');
set(h{4},'Position',FigPosition,'PaperOrientation','landscape','PaperUnits','inches','PaperSize',[8.5*2.1 11*2.1],'PaperPositionMode','auto','Color',Fig_Color)

for i = 1:length(h)-1
    ax_old = get(h{i},'CurrentAxes');
    ax_new = subplot(1,3,i);
    copyobj(allchild(ax_old),ax_new);
    h_lines = allchild(ax_new);
    xmax = zeros(1,5);
    ymax_ln = zeros(1,5);
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
        set(ax_new,'FontSize',FontSize,'XLim',[0 10^(log10(xmax)+0.1)],'YLim',[0 ymax_ax],'XScale','log','XTick',XTick,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth,'Units',Units,'Position',[Axis_Left(i) Axis_Bottom Axis_Width Axis_Height])
    else
        XTick = 0:2:10;
        set(ax_new,'FontSize',FontSize,'XLim',[0 11],'YLim',[0 ymax_ax],'XScale','linear','XTick',XTick,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth,'Units',Units,'Position',[Axis_Left(i) Axis_Bottom Axis_Width Axis_Height])
    end
    title(Title{i})
    xlabel('# of Ambient Dimensions')
    if i == 1
        ylabel('Training Time (sec)')
    end
    hL = legend(ax_new,'Random Forest','R''er F(d)','R''er F(s)','R''er F(s+d)');
    legend(ax_new,'hide')
    get(ax_new,'Position');
end
set(hL,'FontSize',FontSize,'Units',Units,'Position',[Legend_Left Legend_Bottom Legend_Width Legend_Height],'Visible','On','Box',Box_Legend)

fname = '~/LOVEFest/Figures/Fig2_Time';
save_fig(gcf,fname)
