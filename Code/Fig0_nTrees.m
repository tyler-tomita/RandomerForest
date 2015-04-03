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

Colors = linspecer(4,'sequential');
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
Legend_Height = 160;
Legend_Left = 950;
Legend_Bottom = round(Axis_Bottom+Axis_Height/2-Legend_Height/2);
MarkerSize = 14;
Box = 'off';
Box_Legend = 'off';
xmax = [2000 2000 200];

Filename = '~/LOVEFest/Figures/fig/ntrees_to_stabilize.fig';

h_fig_old = openfig(Filename,'invisible');
h_ax_old = flipud(findobj(h_fig_old,'Type','axes'));

h_fig_new = figure('Visible','On');
set(h_fig_new,'Position',FigPosition,'PaperOrientation','landscape','PaperUnits','inches','PaperSize',[8.5*2 11*2],'PaperPositionMode','auto','Color',Fig_Color)

for i = 1:length(h_ax_old)
    h_ax_new = subplot(1,3,i);
    copyobj(allchild(h_ax_old(i)),h_ax_new);
    h_lines = allchild(h_ax_new);
    ymax_ln = zeros(1,length(h_lines));
    for j = 1:length(h_lines)
        set(h_lines(j),'Color',Colors(j,:),'linewidth',LineWidth,'Marker',Marker,'MarkerFaceColor',Colors(j,:),'MarkerEdgeColor',Colors(j,:))
        ymax_ln(j) = max(get(h_lines(j),'YData'));
    end
    ymax_ax = max(ymax_ln);
    ymax_idx = find(ymax_ln==ymax_ax);
    ymax_idx = ymax_idx(1);
    XTick = cat(2,1,xmax(i)/2:xmax(i)/2:xmax(i));
    set(h_ax_new,'XLim',[0 xmax(i)],'YLim',[0 ymax_ax],'XTick',XTick,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth,'Units',Units,'Position',[Axis_Left(i) Axis_Bottom Axis_Width Axis_Height])
    title(Title{i})
    xlabel('Number of Trees')
    ylabel('L hat')
    hL = legend('Random Forest','Dense Randomer Forest','Sparse Randomer Forest','Sparse Randomer Forest w/ Mean Diff');
    legend(h_ax_new,'hide')
    get(h_ax_new,'Position');
end
set(hL,'Units',Units,'Position',[Legend_Left Legend_Bottom Legend_Width Legend_Height],'Visible','On','Box',Box_Legend)

fname = '~/LOVEFest/Figures/Fig0_nTrees';
save_fig(gcf,fname)