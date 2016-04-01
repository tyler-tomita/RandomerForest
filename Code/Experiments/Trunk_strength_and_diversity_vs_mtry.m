close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

Color = [0 1 1;0 1 0;1 0 1;1 0 0;0 0 0];

load Trunk_data
ndims = length(dims);
ntrees = 500;
nmixs = 2:6;
NWorkers = 2;

for i = 4:4

    d = dims(i);
    
    if d <= 5
        mtrys = 1:d;
    elseif d <= 50
        mtrys = ceil(d.^[0 1 2]);
    else
        mtrys = [d.^[0 1] 10*d];
    end
    
    if d >= 6
        nmixs = 2:6;
    else
        nmixs = 2:d;
    end 
    
    for trial = 1:10

        fprintf('trial %d\n',trial)

        for j = 1:length(mtrys)
            
            mtry = mtrys(j);

            fprintf('mtry = %d\n',mtry)
            
            poolobj = gcp('nocreate');
            if isempty(poolobj)
                parpool('local',NWorkers);
            end
            
            tic;
            rerf = rpclassificationforest(ntrees,X{i}(:,:,trial),...
                Y{i}(:,trial),'sparsemethod','sparse',...
                'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            
            Lhat(trial,j) = oobpredict(rerf,X{i}(:,:,trial),Y{i}(:,trial),'last');

            [Yhats,err(j,:)] = oobpredict2(rerf,X{i}(:,:,trial),Y{i}(:,trial));
            dv(trial,j) = diversity(Yhats);
        end
        strength(trial,:) = 1 - mean(err,2)';
    end
end

ax = subplot(3,1,1);
errorbar(mtrys,mean(strength),std(strength)/sqrt(ntrials),'LineWidth',2)
xlabel('mtry')
ylabel('strength')
title(sprintf('Trunk (n = %d, p = %d)',n,d))
ax.XScale = 'log';

ax = subplot(3,1,2);
errorbar(mtrys,mean(dv),std(dv)/sqrt(ntrials),'LineWidth',2)
xlabel('mtry')
ylabel('diversity')
ax.XScale = 'log';

ax = subplot(3,1,3);
errorbar(mtrys,mean(Lhat),std(Lhat)/sqrt(ntrials),'LineWidth',2)
xlabel('mtry')
ylabel('ensemble error rate')
ax.XScale = 'log';

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Trunk_strength_and_diversity_vs_mtry_n' num2str(n) '_d' num2str(d)])