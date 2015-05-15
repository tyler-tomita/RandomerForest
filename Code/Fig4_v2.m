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
FigPosition = [0 140 1300 350];
%left = [50 350 650];
%bottom = [350 50];
%Axis_Left = repmat(left,1,2);
Axis_Left = [210 775];
%Axis_Bottom = cat(2,repmat(bottom(1),1,3),repmat(bottom(2),1,3));
Axis_Bottom = 60;
Axis_Width = 250;
Axis_Height = 250;
Legend_Width = 75;
Legend_Height = 80;
Legend_Left = 475;
Legend_Bottom = round(Axis_Bottom+Axis_Height/2-Legend_Height/2);
MarkerSize = 24;
Box = 'off';
Box_Legend = 'off';
Colorbar_Left = 1100;
Colorbar_Bottom = 30;
Colorbar_Width = 25;
Colorbar_Height = 280;
FontSize = 24;

BasePath = '~/LOVEFest/Figures/fig/';
Filename = {'Fig4_Real_Data_Panel_A_v2.fig','Fig4_Real_Data_Panel_B_v2.fig'};

for i = 1:length(Filename)
    h{i} = openfig(strcat(BasePath,Filename{i}),'invisible');
    grid off
end

h{3} = figure('Visible','On');
set(h{3},'Units','normalized','position',[0 0 1 1]);
set(h{3},'Units','inches');
screenposition = get(h{3},'Position');
set(h{3},...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',screenposition(3:4));

ax_old = get(h{1},'CurrentAxes');
ax_new = subplot(1,2,1);
axis square
copyobj(allchild(ax_old),ax_new);
h_lines = allchild(ax_new);
for j = 1:length(h_lines)
    set(h_lines(j),'Color',Colors(j,:),'linewidth',LineWidth,'Marker',Marker{j},'MarkerSize',MarkerSize,'MarkerFaceColor',Colors(j,:),'MarkerEdgeColor',Colors(j,:))
end
set(ax_new,'FontSize',FontSize,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth)
xlabel('Training Time (sec)')
ylabel('Error Rate')
title('(A) Average Error Rate vs. Time')

ax_old = get(h{2},'CurrentAxes');
ax_new(2) = subplot(1,2,2);
axis square
copyobj(allchild(ax_old),ax_new(2));
h_lines = allchild(ax_new(2));
set(h_lines(1),'Color','k','linewidth',LineWidth)
set(ax_new(2),'FontSize',FontSize,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth)
xlabel('Random Forest')
ylabel('RerF(s+d+r)')
title('(B) Error Rate For Each Dataset')
%hL = legend(ax_new,'Random Forest','Dense Randomer Forest','Sparse Randomer Forest','Sparse Randomer Forest w/ Mean Diff','Robust Sparse Randomer Forest w/ Mean Diff');
legend(ax_new(2),'hide')
get(ax_new(2),'Position');
rgb = map2color(transpose(linspace(1,1000,1000)),'log');
colormap([rgb(:,1) rgb(:,2) rgb(:,3)])
%colorbar
caxis([4 263])
cb = colorbar;
axpos = get(ax_new(end),'Position');
ax_edge = axpos(1) + axpos(3);
cbpos = get(cb,'Position');
set(cb,'position',[ax_edge-.02 axpos(2)+axpos(4)/2-cbpos(4)/2 cbpos(3:4)],'Ticks',[10 25 50 100 200],'FontSize',FontSize)


pos2 = get(ax_new(2),'Position');

pos1 = get(ax_new(1),'Position');
set(ax_new(1),'Position',[pos1(1) pos1(2) pos2(3) pos2(4)])

hL = legend(ax_new(1),'RF','RerF(d)','RerF(s)','RerF(s+d)','RerF(s+d+r)');
axpos = get(ax_new(1),'Position');
ax_edge = axpos(1) + axpos(3);
lpos = get(hL,'Position');
set(hL, 'position',[ax_edge axpos(2)+axpos(4)/2-lpos(4)/2 lpos(3:4)],'Box',Box_Legend,'FontSize',FontSize)
   
pos_l = get(hL,'Position');
left_adjust = pos2(1) - (pos_l(1) + pos_l(3) + .1);

set(ax_new(2),'Position',[pos2(1)-left_adjust/2 pos2(2) pos2(3) pos2(4)])

pos1 = get(ax_new(1),'Position');
set(ax_new(1),'Position',[pos1(1)+left_adjust/2 pos1(2) pos1(3) pos1(4)])

set(hL,'Position',[pos_l(1)+left_adjust/2 pos_l(2) pos_l(3) pos_l(4)])

pos = get(cb,'Position');
set(cb,'Position',[pos(1)-left_adjust/2 pos(2) pos(3) pos(4)])

fname = '~/LOVEFest/Figures/Fig4';
save_fig(gcf,fname)