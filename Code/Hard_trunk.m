%Make Fig1 for LOVEFest (Randomer Forest)

close all
clear
clc

stream = RandStream('mt19937ar','Seed',10);
RandStream.setGlobalStream(stream);

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
Title = {'(A) Trunk'};
Units = 'pixels';
sf = 0.8;
MarkerSize = 24;
Box = 'off';
Box_Legend = 'off';
FontSize = 24;

Level_curve = 1;

BasePath = '~/LOVEFest/Figures/fig/';
Filename = {'Trunk_harder_ooberror_vs_d_n100_ntrees1500_ntrials10.fig'};

%Copy bayes error line to the axes containing Lhat for the different
%classifiers
for i = 1:length(Filename)
    h{i} = openfig(strcat(BasePath,Filename{i}),'invisible');
    grid off
end

h{2} = figure('Visible','On');
set(h{2},'Units','normalized','position',[0 0 1 1]);
set(h{2},'Units','inches');
screenposition = get(h{2},'Position');
set(h{2},...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',screenposition(3:4));

for i = 1:length(h)-1
    ax_old = get(h{i},'CurrentAxes');
    ax_new = gca;
    axis square
    copyobj(allchild(ax_old),ax_new);
    h_lines = allchild(ax_new);
    xmax = zeros(1,length(h_lines));
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
        set(ax_new,'FontSize',FontSize,'XLim',[0 10^(log10(xmax)+0.1)],'YLim',[0 ymax_ax],'XScale','log','XTick',XTick,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth)
    else
        XTick = 0:2:10;
        set(ax_new,'FontSize',FontSize,'XLim',[0 11],'YLim',[0 ymax_ax],'XScale','linear','XTick',XTick,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth)
    end
    title(Title{i})
    xlabel('# Ambient Dim')
    if i == 1
        ylabel('Error Rate')
    end
    a(i) = ax_new;
end
axpos = get(a(end),'Position');
ax_edge = axpos(1) + axpos(3);
hL = legend(a(end),'RF','RerF(d)','RerF(s)','RerF(s+d)','RerF(s+d+r)');
lpos = get(hL,'Position');
set(hL, 'position',[ax_edge axpos(2)+axpos(4)/2-lpos(4)/2 lpos(3:4)],'Box',Box_Legend,'FontSize',FontSize)

pos_l = get(hL(1),'Position');
left_adjust = 1 - (pos_l(1) + pos_l(3));
%overall_width = pos_l(1) + pos_l(3) - pos_a(1);
%if overall_width <= 1
%    left = (1 - overall_width)/2;
%    left_adjust = left - pos_a(1);
%end
for i = 1:length(a)
    pos = get(a(i),'Position');
    set(a(i),'Position',[pos(1)+left_adjust pos(2) pos(3) pos(4)])
end
pos = get(hL(1),'Position');
set(hL(1),'Position',[pos(1)+left_adjust pos(2) pos(3) pos(4)])

fname = '~/LOVEFest/Figures/Trunk_hard';
save_fig(gcf,fname)