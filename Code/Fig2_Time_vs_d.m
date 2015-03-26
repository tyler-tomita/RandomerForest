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
LineWidth = 1.5;
Marker = 'none';
Title = {'Trunk' 'Parity' 'Multimodal'};
Units = 'pixels';
FigPosition = [0 140 1150 325];
Axis_Left = [50 350 650];
Axis_Bottom = 50;
Axis_Width = 250;
Axis_Height = 250;
Legend_Width = 150;
Legend_Height = 80;
Legend_Left = 950;
Legend_Bottom = round(Axis_Bottom+Axis_Height/2-Legend_Height);
MarkerSize = 14;
Box = 'off';
Box_Legend = 'off';

BasePath = '~/LOVEFest/Figures/fig/';
Filename = {'Trunk_time_vs_d_n100_var1_ntrees1000_ntrials10_v2.fig'...
 'Parity_time_vs_d_n100_ntrees1000_ntrials10_v2.fig'...
 'Multimodal_time_vs_d_n100_var1_embed1000_ntrees10_ntrials_v2.fig'};

for i = 1:length(Filename)
    h{i} = openfig(strcat(BasePath,Filename{i}),'invisible');
    grid off
end

h{4} = figure('Visible','On');
set(h{4},'Position',FigPosition,'PaperOrientation','landscape','PaperUnits','inches','PaperSize',[8.5*2 11*2],'PaperPositionMode','auto','Color',Fig_Color)

for i = 1:length(h)-1
    ax_old = get(h{i},'CurrentAxes');
    ax_new = subplot(1,3,i);
    copyobj(allchild(ax_old),ax_new);
    h_lines = allchild(ax_new);
    xmax = zeros(1,5);
    ymax_ln = zeros(1,5);
    for j = 1:4
        set(h_lines(j),'Color',Colors(j,:),'linewidth',LineWidth,'Marker',Marker,'MarkerFaceColor',Colors(j,:),'MarkerEdgeColor',Colors(j,:))
        xmax(j) = max(get(h_lines(j),'XData'));
        ymax_ln(j) = max(get(h_lines(j),'YData'));
    end
    xmax = max(xmax);
    ymax_ax = max(ymax_ln);
    ymax_idx = find(ymax_ln==ymax_ax);
    ymax_ax = ymax_ax + h_lines(ymax_idx).UData(find(h_lines(ymax_idx).YData==ymax_ln(ymax_idx)));
    XTick = logspace(0,log10(xmax),log10(xmax)+1);
    set(ax_new,'XLim',[0 10^(log10(xmax)+0.1)],'YLim',[0 ymax_ax],'XScale','log','XTick',XTick,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth,'Units',Units,'Position',[Axis_Left(i) Axis_Bottom Axis_Width Axis_Height])
    title(Title{i})
    xlabel('Number of Ambient Dimensions')
    ylabel('Training Time (sec)')
    hL = legend(ax_new,'Random Forest','Dense Randomer Forest','Sparse Randomer Forest','Sparse Randomer Forest w/ Mean Diff');
    legend(ax_new,'hide')
    get(ax_new,'Position');
end
set(hL,'Units',Units,'Position',[Legend_Left Legend_Bottom Legend_Width Legend_Height],'Visible','On','Box',Box_Legend)

fname = '~/LOVEFest/Figures/Fig2_Time';
save_fig(gcf,fname)