close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

C = [0 1 1;0 1 0;1 0 1;1 0 0;0 0 0];
Colors.rf = C(1,:);
Colors.rerf = C(2,:);
Colors.rf_rot = C(3,:);
Colors.rerfr = C(4,:);
Colors.frc = C(5,:);
LineWidth = 2;
FontSize = .16;
axWidth = 1.3;
axHeight = 1.3;
% axLeft = [FontSize*4,FontSize*8+axWidth,FontSize*12+axWidth*2,...
%     FontSize*16+axWidth*3,FontSize*4,FontSize*8+axWidth,...
%     FontSize*12+axWidth*2,FontSize*16+axWidth*3];
% axBottom = [FontSize*8+axHeight,FontSize*8+axHeight,FontSize*8+axHeight,...
%     FontSize*8+axHeight,FontSize*4,FontSize*4,FontSize*4,FontSize*4];
axLeft = [FontSize*5,FontSize*6+axWidth,FontSize*7+axWidth*2,...
    FontSize*8+axWidth*3,FontSize*9+axWidth*4,FontSize*5,...
    FontSize*6+axWidth,FontSize*7+axWidth*2,FontSize*8+axWidth*3,...
    FontSize*9+axWidth*4];
axBottom = [FontSize*6+axHeight,FontSize*6+axHeight,FontSize*6+axHeight,...
    FontSize*6+axHeight,FontSize*6+axHeight,FontSize*3,FontSize*3,...
    FontSize*3,FontSize*3 FontSize*3];
% axLeft = [FontSize*4 FontSize*8+axWidth FontSize*4 FontSize*8+axWidth];
% axBottom = [FontSize*6+axHeight FontSize*6+axHeight FontSize*3 FontSize*3];
figWidth = axLeft(end) + axWidth + FontSize*4;
figHeight = axBottom(1) + axHeight + FontSize*3;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

runSims = false;

%% Sparse parity

if runSims
    run_Sparse_parity_transformations
else
    load Sparse_parity_transformations.mat
end

Transformations = fieldnames(mean_err_rf);

for j = 1:length(Transformations)
    Transform = Transformations{j};
    
    % Only keep frc results for nmix = 2
    for i = 1:length(dims)
        d = dims(i);
        if d <= 5
            mtrys = 1:d;
        else
            mtrys = ceil(d.^[0 1/4 1/2 3/4 1]);
        end

        if d >= 6
            nmixs = 2:6;
        else
            nmixs = 2:d;
        end

        for m = 1:length(mtrys)
            mtry = mtrys(m);
            for k = 1:length(nmixs)
                nmix = nmixs(k);
                if nmix ~= 2;
                    mean_err_frc.(Transform)(:,length(nmixs)*(m-1)+k,:) = NaN;
                    sem_err_frc.(Transform)(:,length(nmixs)*(m-1)+k,:) = NaN;
                end
            end
        end
    end
    
    [Lhat.rf,minIdx.rf] = min(mean_err_rf.(Transform)(end,:,:),[],2);
    [Lhat.rerf,minIdx.rerf] = min(mean_err_rerf.(Transform)(end,:,:),[],2);
%     [Lhat.rerfdn,minIdx.rerfdn] = min(mean_err_rerfdn.(Transform)(end,:,:),[],2);
    [Lhat.rf_rot,minIdx.rf_rot] = min(mean_err_rf_rot.(Transform)(end,:,:),[],2);
    [Lhat.rerfr,minIdx.rerfr] = min(mean_err_rerfr.(Transform)(end,:,:),[],2);
    [Lhat.frc,minIdx.frc] = min(mean_err_frc.(Transform)(end,:,:),[],2);

    for i = 1:length(dims)
        sem.rf(i) = sem_rf.(Transform)(end,minIdx.rf(i),i);
        sem.rerf(i) = sem_rerf.(Transform)(end,minIdx.rerf(i),i);
%         sem.rerfdn(i) = sem_rerfdn.(Transform)(end,minIdx.rerfdn(i),i);
        sem.rf_rot(i) = sem_rf_rot.(Transform)(end,minIdx.rf_rot(i),i);
        sem.rerfr(i) = sem_rerfr.(Transform)(end,minIdx.rerfr(i),i);
        sem.frc(i) = sem_frc.(Transform)(end,minIdx.frc(i),i);
    end

    classifiers = fieldnames(Lhat);
    
    ax = subplot(2,5,j);
    
    for i = 1:length(classifiers)
        cl = classifiers{i};
        h = errorbar(dims,Lhat.(cl)(:)',sem.(cl),'LineWidth',LineWidth,'Color',Colors.(cl));
        hold on
    end
    
    title(['(' char('A'+j-1) ')'],'Units','normalized','Position',[0.025 0.975],'HorizontalAlignment','left','VerticalAlignment','top')
    text(0.5,1,Transform,'FontSize',14,'FontWeight','bold','Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom')
%     title(['(' char('A'+j-1) ') ' Transform])
    if j == 1
        xlabel('p')
        ylabel({'\bf{Sparse Parity}';'\rm{Error Rate}'})
        ax.XTick = [5 10 25 50 100];
        ax.XTickLabel = {'5';'10';'25';'50';'100'};
    else
        ax.XTick = [];
        ax.YTickLabel = {};
    end
    ax.LineWidth = LineWidth;
    ax.FontUnits = 'inches';
    ax.FontSize = FontSize;
    ax.Units = 'inches';
    ax.Position = [axLeft(j) axBottom(j) axWidth axHeight];
    ax.Box = 'off';
    ax.XLim = [0 105];
    ax.XScale = 'log';
    ax.YLim = [0 .6];
    ax.YGrid = 'on';
    ax.YMinorGrid = 'on';
    ax.GridLineStyle = ':';
end


clear Lhat sem minIdx

%% Trunk

if runSims
    run_Trunk_transformations
else
    load Trunk_transformations2.mat
end

Transformations = fieldnames(mean_err_rf);

for j = 1:length(Transformations)
    Transform = Transformations{j};
    
        % Only keep frc results for nmix = 2
    for i = 1:length(dims)
        d = dims(i);
        if d <= 5
            mtrys = 1:d;
        else
            mtrys = ceil(d.^[0 1/4 1/2 3/4 1]);
        end

        if d >= 6
            nmixs = 2:6;
        else
            nmixs = 2:d;
        end

        for m = 1:length(mtrys)
            mtry = mtrys(m);
            for k = 1:length(nmixs)
                nmix = nmixs(k);
                if nmix ~= 2;
                    mean_err_frc.(Transform)(:,length(nmixs)*(m-1)+k,:) = NaN;
                    sem_err_frc.(Transform)(:,length(nmixs)*(m-1)+k,:) = NaN;
                end
            end
        end
    end
    
    [Lhat.rf,minIdx.rf] = min(mean_err_rf.(Transform)(end,:,:),[],2);
    [Lhat.rerf,minIdx.rerf] = min(mean_err_rerf.(Transform)(end,:,:),[],2);
%     [Lhat.rerfdn,minIdx.rerfdn] = min(mean_err_rerfdn.(Transform)(end,:,:),[],2);
    [Lhat.rf_rot,minIdx.rf_rot] = min(mean_err_rf_rot.(Transform)(end,:,:),[],2);
    [Lhat.rerfr,minIdx.rerfr] = min(mean_err_rerfr.(Transform)(end,:,:),[],2);
    [Lhat.frc,minIdx.frc] = min(mean_err_frc.(Transform)(end,:,:),[],2);

    for i = 1:length(dims)
        sem.rf(i) = sem_rf.(Transform)(end,minIdx.rf(i),i);
        sem.rerf(i) = sem_rerf.(Transform)(end,minIdx.rerf(i),i);
%         sem.rerfdn(i) = sem_rerfdn.(Transform)(end,minIdx.rerfdn(i),i);
        sem.rf_rot(i) = sem_rf_rot.(Transform)(end,minIdx.rf_rot(i),i);
        sem.rerfr(i) = sem_rerfr.(Transform)(end,minIdx.rerfr(i),i);
        sem.frc(i) = sem_frc.(Transform)(end,minIdx.frc(i),i);
    end

    classifiers = fieldnames(Lhat);
    
    ax = subplot(2,5,j+5);
    
    for i = 1:length(classifiers)
        cl = classifiers{i};
        h = errorbar(dims,Lhat.(cl)(:)',sem.(cl),'LineWidth',LineWidth,'Color',Colors.(cl));
        hold on
    end
    
    title(['(' char('A'+j+4) ')'],'Units','normalized','Position',[0.025 0.975],'HorizontalAlignment','left','VerticalAlignment','top')
%     text(0.5,1,Transform,'FontSize',14,'FontWeight','bold','Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom')
%     title(['(' char('A'+j+4) ') ' Transform])
    ax.LineWidth = LineWidth;
    ax.FontUnits = 'inches';
    ax.FontSize = FontSize;
    ax.Units = 'inches';
    ax.Position = [axLeft(j+5) axBottom(j+5) axWidth axHeight];
    ax.Box = 'off';
    ax.XLim = [1 600];
    ax.YLim = [0.02 .17];
    ax.YGrid = 'on';
    ax.YMinorGrid = 'on';
    ax.GridLineStyle = ':';
    ax.XScale = 'log';
    if j+5 == 6
        xlabel('p')
        ylabel({'\bf{Trunk}';'\rm{Error Rate}'})
        ax.XTick = [logspace(0,2,3) 500];
        ax.XTickLabel = {'1';'10';'100';'500'};
    else
        ax.XTick = [];
        ax.YTickLabel = {};
    end
end

% l = legend('RF','RerF','RotRF');
% l.Location = 'southeast';
% l.Box = 'off';
% l.FontSize = 10;

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Fig3_transformations2'])