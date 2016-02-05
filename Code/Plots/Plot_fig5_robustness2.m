close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

C = [0 1 1;0 1 0;1 0 1;1 0 0;0 0 0;1 .5 0];
LineWidth = 2;
FontSize = .16;
axWidth = 1.3;
axHeight = 1.3;
axLeft = repmat(FontSize*5,5,1);
axBottom = [FontSize*15+axHeight*4;FontSize*12+axHeight*3;FontSize*9+axHeight*2;FontSize*6+axHeight;FontSize*3];
figWidth = axLeft(end) + axWidth + FontSize*12;
figHeight = axBottom(1) + axHeight + FontSize*3;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

Titles = {'RF' 'RerF' 'RerF(d)' 'RerF(d+r)' 'Rotation RF'};

runSims = false;

if runSims
    run_Trunk_transformations
else
    load Trunk_transformations.mat
    load Trunk_transformations_rankRerFd.mat
end

Transformations = fieldnames(mean_err_rf);

for j = 1:length(Transformations)
    Transform = Transformations{j};
    
    [Lhat.rf,minIdx.rf] = min(mean_err_rf.(Transform)(end,:,:),[],2);
    [Lhat.rerf,minIdx.rerf] = min(mean_err_rerf.(Transform)(end,:,:),[],2);
    [Lhat.rerfdn,minIdx.rerfdn] = min(mean_err_rerfdn.(Transform)(end,:,:),[],2);
    [Lhat.rerfdnr,minIdx.rerfdnr] = min(mean_err_rerfdnr.(Transform)(end,:,:),[],2);
    [Lhat.rf_rot,minIdx.rf_rot] = min(mean_err_rf_rot.(Transform)(end,:,:),[],2);

    for i = 1:length(dims)
        sem.rf(i) = sem_rf.(Transform)(end,minIdx.rf(i),i);
        sem.rerf(i) = sem_rerf.(Transform)(end,minIdx.rerf(i),i);
        sem.rerfdn(i) = sem_rerfdn.(Transform)(end,minIdx.rerfdn(i),i);
        sem.rerfdnr(i) = sem_rerfdnr.(Transform)(end,minIdx.rerfdnr(i),i);
        sem.rf_rot(i) = sem_rf_rot.(Transform)(end,minIdx.rf_rot(i),i);
    end

    classifiers = fieldnames(Lhat);
    
    for i = 1:length(classifiers)
        cl = classifiers{i};
        if j == 1
            ax(i) = subplot(5,1,i);
            hold on
        else
            axes(ax(i));
        end
        h = errorbar(dims,Lhat.(cl)(:)',sem.(cl),'LineWidth',LineWidth,'Color',C(j,:));
    end  
end

for i = 1:length(classifiers)
    axes(ax(i));
    title(['(' char('A'+i-1) ')'],'Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
    text(-0.55,0.5,Titles{i},'FontSize',14,'FontWeight','bold',...
        'Units','normalized','HorizontalAlignment','center',...
        'VerticalAlignment','middle','Rotation',90)
    if i == 1
        xlabel('d')
        ylabel('Error Rate')
    end
    ax(i).LineWidth = LineWidth;
    ax(i).FontUnits = 'inches';
    ax(i).FontSize = FontSize;
    ax(i).Units = 'inches';
    ax(i).XLim = [1 550];
    ax(i).XScale = 'log';
    ax(i).XTick = [logspace(0,2,3) 500];
    ax(i).YLim = [0 0.175];
    if i > 1
        ax(i).XTickLabel = {};
        ax(i).YTickLabel = {};
    else
        ax(i).XTickLabel = {'1' '10' '100' '500'};
    end
    ax(i).Position = [axLeft(i) axBottom(i) axWidth axHeight];
    ax(i).Box = 'off';
end

l = legend(Transformations);
l.Location = 'northeast';
l.Box = 'off';
l.FontSize = 10;

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Fig5_robustness2'])