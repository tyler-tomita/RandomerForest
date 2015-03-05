%Make Fig1 for LOVEFest (Randomer Forest)

%Outline:
%Loop over the three existing figures
%Format the figure properties to make it aesthetically pleasing
%Copy to 1x3 subplot

close all
clear
clc

Colors = linspecer(5,'sequential');
Fig_Color = [1 1 1];
LineWidth = 1.5;
Marker = 'none';
Title = {'Trunk' 'Parity' 'Multimodal'};
Units = 'pixels';
FigPosition = [0 140 1150 650];
Axis_Left = [50 350 650 50 350 650];
Axis_Bottom = [350 350 350 50 50 50];
Axis_Width = 250;
Axis_Height = 250;
Legend_Left = 950;
Legend_Bottom = [475 175];
Legend_Width = 150;
Legend_Height = 100;
MarkerSize = 14;
Box = 'off';

BasePath = '~/LOVEFest/Figures/fig/';
Filename = {'trunk_time_vs_d_n100_var1_embed100_ntrees1000_ntrials10.fig'...
 'Parity_time_vs_d_n100_var1_embed100_ntrees1000_ntrials10.fig'...
 'Multimodal_time_vs_d_n100_var1_embed1000_ntrees10_ntrials.fig'};


%Copy bayes error line to the axes containing Lhat for the different
%classifiers
for i = 1:length(Filename)
    h{i} = openfig(strcat(BasePath,Filename{i}),'invisible');
    grid off
end

h{4} = figure('Visible','On');
set(h{4},'Position',FigPosition,'Color',Fig_Color)

for i = 1:length(h)-1
    ax_old = get(h{i},'CurrentAxes');
    ax_new = subplot(1,3,i);
    copyobj(allchild(ax_old),ax_new);
    h_lines = allchild(ax_new);
    xmax = zeros(1,4);
    ymax = zeros(1,4);
    for j = 1:4
        set(h_lines(j),'Color',Colors(j,:),'linewidth',LineWidth,'Marker',Marker,'MarkerFaceColor',Colors(j,:),'MarkerEdgeColor',Colors(j,:))
        xmax(j) = max(get(h_lines(j),'XData'));
        ymax(j) = max(get(h_lines(j),'YData'));
    end
    xmax = max(xmax);
    ymax = max(ymax);
    XTick = logspace(0,log10(xmax),log10(xmax)+1);
    set(ax_new,'XLim',[0 xmax],'YLim',[0 ymax],'XScale','log','XTick',XTick,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth,'Units',Units,'Position',[Axis_Left(i) Axis_Bottom(i) Axis_Width Axis_Height])
    title(Title{i})
    xlabel('Number of Ambient Dimensions')
    ylabel('Time')
    hL = legend(ax_new,'Random Forest','Dense Randomer Forest','Sparse Randomer Forest','Sparse Randomer Forest w/ Mean Diff');
    legend(ax_new,'hide')
    get(ax_new,'Position');
end
set(hL,'Units',Units,'Position',[Legend_Left Legend_Bottom(1) Legend_Width Legend_Height],'Visible','On','Box','Off')
fname = 'Fig0_Time';
save_fig(gcf,fname)