close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

LineWidth = 2;
FontSize = .2;
axWidth = 1.5;
axHeight = 1.5;
axLeft = [FontSize*4,FontSize*8+axWidth,FontSize*12+axWidth*2,...
    FontSize*16+axWidth*3,FontSize*4,FontSize*8+axWidth,...
    FontSize*12+axWidth*2,FontSize*16+axWidth*3];
axBottom = [FontSize*8+axHeight,FontSize*8+axHeight,FontSize*8+axHeight,...
    FontSize*8+axHeight,FontSize*4,FontSize*4,FontSize*4,FontSize*4];
figWidth = axLeft(end) + FontSize*4;
figHeight = axBottom(1) + axHeight + FontSize*4;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];

Titles = {'(A) RF' '(B) RerF' '(C) RerFd' '(D) Rotation RF',...
    '(E) RF' '(F) RerF' '(G) RerFd' '(H) Rotation RF'};

runSims = false;

if runSims
    run_Sparse_parity_transformations
else
    load Sparse_parity_transformations.mat
end

Transformations = fieldnames(mean_err_rf);

for j = 1:length(Transformations)
    Transform = Transformations{j};
    
    [Lhat.rf,minIdx.rf] = min(mean_err_rf.(Transform)(end,:,:),[],2);
    [Lhat.rerf,minIdx.rerf] = min(mean_err_rerf.(Transform)(end,:,:),[],2);
    [Lhat.rerfdn,minIdx.rerfdn] = min(mean_err_rerfdn.(Transform)(end,:,:),[],2);
    [Lhat.rf_rot,minIdx.rf_rot] = min(mean_err_rf_rot.(Transform)(end,:,:),[],2);

    for i = 1:length(dims)
        sem.rf(i) = sem_rf.(Transform)(end,minIdx.rf(i),i);
        sem.rerf(i) = sem_rerf.(Transform)(end,minIdx.rerf(i),i);
        sem.rerfdn(i) = sem_rerfdn.(Transform)(end,minIdx.rerfdn(i),i);
        sem.rf_rot(i) = sem_rf_rot.(Transform)(end,minIdx.rf_rot(i),i);
    end

    classifiers = fieldnames(Lhat);
    
    for i = 1:length(classifiers)
        cl = classifiers{i};
        if j == 1
            ax(i) = subplot(2,4,i);
            hold on
        else
            axes(ax(i));
        end
        h = errorbar(dims,Lhat.(cl)(:)',sem.(cl),'LineWidth',LineWidth);
    end  
end

for i = 1:length(classifiers)
    axes(ax(i));
    title(Titles{i})
    xlabel('d')
    if i == 1
        ylabel('Error Rate')
    end
    ax(i).LineWidth = LineWidth;
    ax(i).FontUnits = 'inches';
    ax(i).FontSize = FontSize;
    ax(i).Units = 'inches';
    ax(i).Position = [axLeft(i) axBottom(i) axWidth axHeight];
    ax(i).Box = 'off';
end

clear Lhat sem minIdx

if runSims
    run_Trunk_transformations
else
    load Trunk_transformations.mat
end

Transformations = fieldnames(mean_err_rf);

for j = 1:length(Transformations)
    Transform = Transformations{j};
    
    [Lhat.rf,minIdx.rf] = min(mean_err_rf.(Transform)(end,:,:),[],2);
    [Lhat.rerf,minIdx.rerf] = min(mean_err_rerf.(Transform)(end,:,:),[],2);
    [Lhat.rerfdn,minIdx.rerfdn] = min(mean_err_rerfdn.(Transform)(end,:,:),[],2);
    [Lhat.rf_rot,minIdx.rf_rot] = min(mean_err_rf_rot.(Transform)(end,:,:),[],2);

    for i = 1:length(dims)
        sem.rf(i) = sem_rf.(Transform)(end,minIdx.rf(i),i);
        sem.rerf(i) = sem_rerf.(Transform)(end,minIdx.rerf(i),i);
        sem.rerfdn(i) = sem_rerfdn.(Transform)(end,minIdx.rerfdn(i),i);
        sem.rf_rot(i) = sem_rf_rot.(Transform)(end,minIdx.rf_rot(i),i);
    end

    classifiers = fieldnames(Lhat);
    
    for i = 1:length(classifiers)
        cl = classifiers{i};
        if j == 1
            ax(i+4) = subplot(2,4,i+4);
            hold on
        else
            axes(ax(i+4));
        end
        h = errorbar(dims,Lhat.(cl)(:)',sem.(cl),'LineWidth',LineWidth);
    end  
end

for i = 1:length(classifiers)
    axes(ax(i+4));
    title(Titles{i+4})
    xlabel('d')
    if i == 1
        ylabel('Error Rate')
    end
    ax(i+4).LineWidth = LineWidth;
    ax(i+4).FontUnits = 'inches';
    ax(i+4).FontSize = FontSize;
    ax(i+4).Units = 'inches';
    ax(i+4).Position = [axLeft(i+4) axBottom(i+4) axWidth axHeight];
    ax(i+4).Box = 'off';
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Fig3_transformations'])