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
Colors = linspecer(2,'sequential');
Fig_Color = [1 1 1];
LineWidth = 3;
Marker = 'none';
Title = {'(A) RF' '(B) RerF(d)' '(C) RerF(d+r)' '(D) RF' '(E) RerF(d)' '(F) RerF(d+r)'};
Units = 'pixels';
%FigPosition = [0 140 1150 650];
FigPosition = [0 140 1300 350];
%left = [50 350 650];
%bottom = [350 50];
%Axis_Left = repmat(left,1,2);
Axis_Left = [85 410 735];
%Axis_Bottom = cat(2,repmat(bottom(1),1,3),repmat(bottom(2),1,3));
Axis_Bottom = 63;
Axis_Width = 250;
Axis_Height = 250;
Legend_Width = 150;
Legend_Height = 160;
Legend_Left = 1075;
%Legend_Bottom = round(bottom(1)+Axis_Height/2-Legend_Height(1)/2);
Legend_Bottom = round(Axis_Bottom+Axis_Height/2-Legend_Height/2);
MarkerSize = 24;
Box = 'off';
Box_Legend = 'off';
FontSize = 24;
XLim = [-8 8];
YLim = XLim;

Filename = '~/LOVEFest/Figures/fig/Cigars_rough.fig';

h_fig_old = openfig(Filename,'invisible');
ax_old = flipud(findobj(h_fig_old,'Type','axes'));

h_fig_new = figure('Visible','On');
set(h_fig_new,'Units','normalized','position',[0 0 1 1]);
set(h_fig_new,'Units','inches');
screenposition = get(h_fig_new,'Position');
set(h_fig_new,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',screenposition(3:4));

for i = 1:length(ax_old)
    %h_ax_new = subplot(2,3,i);
    ax_new(i) = subplot(1,3,i);
    axis square
    copyobj(allchild(ax_old(i)),ax_new(i));
    h_lines = allchild(ax_new(i));
    if i == 1
        hh = h_lines;
    end
    for j = 1:length(h_lines)
        if j ~= length(h_lines) && j ~= length(h_lines)-1
            set(h_lines(j),'Color','k','linewidth',2,'Marker','none')
        elseif j == length(h_lines)
            set(h_lines(j),'linewidth',LineWidth)
        else
            set(h_lines(j),'linewidth',LineWidth)
        end
    end
    set(ax_new(i),'FontSize',FontSize,'Box',Box,'LineWidth',LineWidth)
    if i < 3
        xlim(XLim)
        ylim(YLim)
    else
        xlim([-500 1500])
        ylim([-500 1500])
    end
    title(Title{i})
    xlabel('Ambient Dimension 1')
    ylabel('Ambient Dimension 2')
end

axpos = get(ax_new(end),'Position');
ax_edge = axpos(1) + axpos(3);

fname = '~/LOVEFest/Figures/Cigars';
save_fig(gcf,fname)
