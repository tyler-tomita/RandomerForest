%% Plot posterior heat maps

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

Colors.rf = 'c';
Colors.rfr = 'c';
Colors.rerf = 'g';
Colors.rerfr = 'g';
Colors.rerf2 = 'b';
Colors.rerf2r = 'b';
Colors.frc = 'm';
Colors.frcr = 'm';
LineStyles.rf = '-';
LineStyles.rfr = ':';
LineStyles.rerf = '-';
LineStyles.rerfr = ':';
LineStyles.rerf2 = '-';
LineStyles.rerf2r = ':';
LineStyles.frc = '-';
LineStyles.frcr = ':';

LineWidth = 2;
FontSize = .15;
axWidth = .75;
axHeight = .75;
cbWidth = .1;
cbHeight = axHeight;
axLeft = repmat([FontSize*3,FontSize*5+axWidth,FontSize*7+2*axWidth],1,8);
% axLeft = [FontSize*3,FontSize*5+axWidth,FontSize*7+2*axWidth,FontSize*3,...
%     FontSize*5+axWidth,FontSize*7+2*axWidth,FontSize*3,...
%     FontSize*5+axWidth,FontSize*7+2*axWidth,FontSize*3,...
%     FontSize*5+axWidth,FontSize*7+2*axWidth,FontSize*3,...
%     FontSize*5+axWidth,FontSize*7+2*axWidth];
axBottom = [(FontSize*13.5+axHeight*8)*ones(1,3),(FontSize*12+axHeight*7)*ones(1,3),(FontSize*10.5+axHeight*6)*ones(1,3),(FontSize*9+axHeight*5)*ones(1,3),...
    (FontSize*7.5+axHeight*4)*ones(1,3),(FontSize*6+axHeight*3)*ones(1,3),...
    (FontSize*4.5+axHeight*2)*ones(1,3),(FontSize*3+axHeight)*ones(1,3)];
axHeight2 = axHeight;
axWidth2 = axLeft(end) - axLeft(1);
axLeft2 = axLeft(1) + axWidth(1)/2;
axBottom2 = FontSize*1.5;
cbLeft = axLeft + axWidth + FontSize/2;
cbBottom = axBottom;
legWidth = axWidth;
legHeight = axHeight2;
legLeft = axLeft2 + axWidth2 + FontSize/2;
legBottom = axBottom2;
figWidth = cbLeft(end) + cbWidth + FontSize*3;
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
Posteriors{3}.truth.Affine = truth.posteriors;
Posteriors{3}.truth.Outlier = truth.posteriors;
Posteriors{3} = orderfields(Posteriors{3},[length(fieldnames(Posteriors{3})),1:length(fieldnames(Posteriors{3}))-1]);

Classifiers = fieldnames(Posteriors{3});

Classifiers(~ismember(Classifiers,{'truth' 'rf','rerf','rerfr','rerf2','rerf2r','frc','frcr'})) = [];

ClassifierNames = {'Posterior' 'RF' 'RerF' 'RerF(r)' 'RerF2','RerF2(r)' 'F-RC' 'Frank'};

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
        ph{(i-1)*3+j} = posterior_map(Xpost,Ypost,mean(Posteriors{3}.(Classifiers{i}).(Transformations{j}),3),false);
%         title(['(' char('A'+(i-1)*3+j-1) ')'],'FontSize',16,'Units','normalized','Position',[-0.02 1],...
%             'HorizontalAlignment','right','VerticalAlignment','top')
        if i==1
            if strcmp(Transformations{j},'Untransformed')
                text(0.5,1.05,'Raw','FontSize',14,'FontWeight','bold','Units',...
                    'normalized','HorizontalAlignment','center','VerticalAlignment'...
                    ,'bottom')
            elseif strcmp(Transformations{j},'Outlier')
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
        %Compute hellinger distance
            ntrials = size(Posteriors{3}.(Classifiers{i}).(Transformations{j}),3);
            for trial = 1:ntrials
                hd.(Classifiers{i})(trial,j) = mean(1/sqrt(2)*sqrt(sum((sqrt(Posteriors{3}.(Classifiers{i}).(Transformations{j})(:,:,trial)) - sqrt(Posteriors{3}.truth.(Transformations{j}))).^2,2)));
            end
%             text(0.5,1.05,sprintf('HD = %0.3f',hd.(Classifiers{i})(trial,j)),'FontSize',16,'Units',...
%                 'normalized','HorizontalAlignment','center','VerticalAlignment'...
%                 ,'bottom')
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
            cb.FontSize = 14;
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

ax(19) = axes;

Classifiers = fieldnames(hd);

for c = 1:length(Classifiers)
errorbar(mean(hd.(Classifiers{c})),std(hd.(Classifiers{c}))/sqrt(ntrials),...
    'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}),...
    'LineStyle',LineStyles.(Classifiers{c}))
hold on
end

ylabel({'Hellinger';'Distance'})
% title('(S)','FontSize',16,'Units','normalized','Position',[-0.02 1],...
%     'HorizontalAlignment','right','VerticalAlignment','top')
ax(19).LineWidth = LineWidth;
ax(19).FontUnits = 'inches';
ax(19).FontSize = FontSize;
ax(19).Units = 'inches';
ax(19).Position = [axLeft2 axBottom2 axWidth2 axHeight2];
ax(19).XLim = [0.9,3.1];
ax(19).YLim = [0.25,0.54];
ax(19).XTick = 1:3;
ax(19).XTickLabel = {'Raw','Affine','Corrupted'};
ax(19).Box = 'off';

[lh,objh] = legend('RF','RerF','RerF(r)','RerF2','RerF2(r)','F-RC','Frank');
lh.Units = 'inches';
lh.Position = [legLeft legBottom legWidth legHeight];
lh.Box = 'off';

for i = 8:length(objh)
    objh(i).Children.Children(2).XData = [(objh(i).Children.Children(2).XData(2)-objh(i).Children.Children(2).XData(1))*0.5+objh(i).Children.Children(2).XData(1),objh(i).Children.Children(2).XData(2)];
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Fig1_posteriors'])