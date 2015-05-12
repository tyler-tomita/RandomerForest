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

%Colors = linspecer(3,'sequential');
Colors = repmat(linspecer(5,'sequential'),2,1);
Fig_Color = [1 1 1];
LineWidth = 3;
Marker = {'none' 'none' 'none' 'none' 'none' '.' '.' '.' '.' '.'};
Units = 'pixels';
%FigPosition = [0 140 1150 650];
FigPosition = [0 140 975 350];
%left = [50 350 650];
%bottom = [350 50];
%Axis_Left = repmat(left,1,2);
Axis_Left = [75 600];
%Axis_Bottom = cat(2,repmat(bottom(1),1,3),repmat(bottom(2),1,3));
Axis_Bottom = 50;
Axis_Width = 250;
Axis_Height = 250;
Legend_Width = 150;
Legend_Height = 80;
Legend_Left = 350;
Legend_Bottom = round(Axis_Bottom+Axis_Height/2-Legend_Height/2);
MarkerSize = 14;
Box = 'off';
Box_Legend = 'off';
Colorbar_Left = 875;
Colorbar_Bottom = 50;
Colorbar_Width = 25;
Colorbar_Height = 250;
FontSize = 20;

BasePath = '~/LOVEFest/Figures/fig/';
Filename = {'Fig4_Real_Data_Panel_A_v2.fig','Fig4_Real_Data_Panel_B_v2.fig'};

for i = 1:length(Filename)
    h{i} = openfig(strcat(BasePath,Filename{i}),'invisible');
    grid off
end

h{3} = figure('Visible','On');
set(h{3},'Position',FigPosition,'PaperOrientation','landscape','PaperUnits','inches','PaperSize',[8.5*2.1 11*2.1],'PaperPositionMode','auto','Color',Fig_Color)

ax_old = get(h{1},'CurrentAxes');
ax_new = subplot(1,2,1);
copyobj(allchild(ax_old),ax_new);
h_lines = allchild(ax_new);
for j = 1:length(h_lines)
    set(h_lines(j),'Color',Colors(j,:),'linewidth',LineWidth,'Marker',Marker{j},'MarkerSize',MarkerSize,'MarkerFaceColor',Colors(j,:),'MarkerEdgeColor',Colors(j,:))
end
set(ax_new,'FontSize',FontSize,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth,'Units',Units,'Position',[Axis_Left(1) Axis_Bottom Axis_Width Axis_Height])
xlabel('Training Time (sec)')
ylabel('Lhat')
title('(A) Average Lhat vs. Time')
hL = legend('Random Forest','R''er F(d)','R''er F(s)','R''er F(s+d)','R''er F(s+d+r)');
set(hL,'Units',Units,'Position',[Legend_Left Legend_Bottom Legend_Width Legend_Height],'Visible','On','Box',Box_Legend)
get(ax_new,'Position');

ax_old = get(h{2},'CurrentAxes');
ax_new = subplot(1,2,2);
copyobj(allchild(ax_old),ax_new);
h_lines = allchild(ax_new);
set(h_lines(1),'Color','k','linewidth',LineWidth)
set(ax_new,'FontSize',FontSize,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth,'Units',Units,'Position',[Axis_Left(2) Axis_Bottom Axis_Width Axis_Height])
xlabel('Random Forest')
ylabel('R''er F(s+d+r)')
title('(B) Lhat For Each Datasets')
%hL = legend(ax_new,'Random Forest','Dense Randomer Forest','Sparse Randomer Forest','Sparse Randomer Forest w/ Mean Diff','Robust Sparse Randomer Forest w/ Mean Diff');
legend(ax_new,'hide')
get(ax_new,'Position');
rgb = map2color(transpose(linspace(1,1000,1000)),'log');
colormap([rgb(:,1) rgb(:,2) rgb(:,3)])
colorbar
caxis([4 263])
h = colorbar('Units',Units,'Position',[Colorbar_Left Colorbar_Bottom Colorbar_Width Colorbar_Height]);
   

fname = '~/LOVEFest/Figures/Fig4';
save_fig(gcf,fname)