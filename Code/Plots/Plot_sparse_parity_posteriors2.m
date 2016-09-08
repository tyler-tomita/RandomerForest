%% Plot posterior heat maps

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);
LineWidth = 2;
FontSize = .2;
axWidth = 1.3;
axHeight = 1.3;
cbWidth = .1;
cbHeight = axHeight;
axLeft = [FontSize*3,FontSize*5+axWidth,FontSize*7+2*axWidth,FontSize*3,...
    FontSize*5+axWidth,FontSize*7+2*axWidth,FontSize*3,...
    FontSize*5+axWidth,FontSize*7+2*axWidth,FontSize*3,...
    FontSize*5+axWidth,FontSize*7+2*axWidth,FontSize*3,...
    FontSize*5+axWidth,FontSize*7+2*axWidth];
axBottom = [(FontSize*7.5+axHeight*4)*ones(1,3),(FontSize*6+axHeight*3)*ones(1,3),...
    (FontSize*4.5+axHeight*2)*ones(1,3),(FontSize*3+axHeight)*ones(1,3),...
    FontSize*1.5*ones(1,3)];
cbLeft = axLeft + axWidth + FontSize/2;
cbBottom = axBottom;
figWidth = cbLeft(end) + cbWidth + FontSize*2;
figHeight = axBottom(1) + axHeight + FontSize*2;

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

runSims = false;

if runSims
    run_sparse_parity_posteriors
else
    load Sparse_parity_true_posteriors
    load Sparse_parity_uniform_transformations_posteriors
end

Posteriors = Phats;
clear Phats

Posteriors{3}.truth.Untransformed = truth.posteriors;
Posteriors{3}.truth.Scaled = truth.posteriors;
Posteriors{3}.truth.Outlier = truth.posteriors;
Posteriors{3} = orderfields(Posteriors{3},[length(fieldnames(Posteriors{3})),1:length(fieldnames(Posteriors{3}))-1]);

Classifiers = fieldnames(Posteriors{3});
Classifiers(strcmp(Classifiers,'rr_rfr')) = [];

ClassifierNames = {'Posterior' 'RF' 'F-RC' 'F-RC(r)' 'RR-RF'};

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

for i = 1:length(Classifiers)
    Transformations = fieldnames(Posteriors{3}.(Classifiers{i}));
    for j = 1:length(Transformations)
        ax((i-1)*3+j) = axes;
        ph{(i-1)*3+j} = posterior_map(Xpost.(Transformations{j}),Ypost.(Transformations{j}),mean(Posteriors{3}.(Classifiers{i}).(Transformations{j}),3),false);
        title(['(' char('A'+(i-1)*3+j-1) ')'],'FontSize',14,'Units','normalized','Position',[-0.02 1],...
            'HorizontalAlignment','right','VerticalAlignment','top')
        if i==1
            if strcmp(Transformations{j},'Outlier')
                text(0.5,1.05,'Corrupted','FontSize',14,'FontWeight','bold','Units',...
                    'normalized','HorizontalAlignment','center','VerticalAlignment'...
                    ,'bottom')
            else
                text(0.5,1.05,Transformations{j},'FontSize',14,'FontWeight','bold','Units',...
                    'normalized','HorizontalAlignment','center','VerticalAlignment'...
                    ,'bottom')
            end
            if j==1
                xlabel('X_1')
                ylabel({['\bf{',ClassifierNames{i},'}'];'\rm{X_2}'})
            end
        else
            if j==1
                ylabel(['\bf{',ClassifierNames{i},'}'])
            end
        end
        ax((i-1)*3+j).LineWidth = LineWidth;
        ax((i-1)*3+j).FontUnits = 'inches';
        ax((i-1)*3+j).FontSize = FontSize;
        ax((i-1)*3+j).Units = 'inches';
        ax((i-1)*3+j).Position = [axLeft((i-1)*3+j) axBottom((i-1)*3+j) axWidth axHeight];
        ax((i-1)*3+j).XTick = [];
        ax((i-1)*3+j).YTick = ax((i-1)*3+j).XTick;
        if i==1
            colormap(ax((i-1)*3+j),'jet')
        else
            colormap(ax((i-1)*3+j),'parula')
        end
        
        if i==1 && j==length(Transformations) || i==length(Classifiers) && j==length(Transformations)
            cb = colorbar;
            cb.Units = 'inches';
            cb.Position = [cbLeft((i-1)*3+j) cbBottom((i-1)*3+j) cbWidth cbHeight];
            cb.Box = 'off';
        end
    end
end

cdata = [];
for i = 4:length(ph)
    cdata = [cdata;ph{i}.CData(:)];
end
cmin = min(cdata);
cmax = max(cdata);

for i = 4:length(ax)
    caxis(ax(i),[cmin cmax])
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Fig1_posteriors'])