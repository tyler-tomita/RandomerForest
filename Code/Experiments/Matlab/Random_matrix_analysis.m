%% Plot number proportion of columns that are all zero and proportion of columns that are nonunique as a function of p and d

close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(123);

ps = [2,5,10,25,50,100,250,500,1000];
slope = zeros(1,length(ps));

for i = 1:length(ps)
    p = ps(i);
    
    fprintf('p = %d\n',p)
    
    if p <= 50
        ntrials = 200;
    else
        ntrials = 50;
    end
    
    if p >= 5 && p <= 500
        ds = ceil(p.^[0 0.5 1 1.5 2]);
    elseif p > 500
        ds = ceil(p.^[0 0.5 1 1.5]);
    else
        ds = [1:p ceil(p.^[1.5 2])];
    end
    nz = zeros(ntrials,length(ds));
    ndup = zeros(ntrials,length(ds));
    dprime = zeros(ntrials,length(ds));
    
    for j = 1:length(ds)
        d = ds(j);
        for trial = 1:ntrials
            fprintf('trial %d\n',trial)
            M = randmat(p,d,'sparse2',1/p);
            dnz = size(M(:,any(M)),2);
            dprime(trial,j) = size(unique(M(:,any(M))','rows')',2);
            nz(trial,j) = (d - dnz)/d;
            ndup(trial,j) = (dnz - dprime(trial,j))/d;
        end
    end
    figure(i)
    ax = subplot(1,2,1);
    hold on
    errorbar(ds,mean(nz),std(nz)/sqrt(ntrials),'LineWidth',2)
    errorbar(ds,mean(ndup),std(ndup)/sqrt(ntrials),'LineWidth',2)
    xlabel('d')
    ylabel('proportion')
    title(sprintf('p = %d',p))
    legend('all-zero columns','duplicate columns','Location','southeast')
    ax.XScale = 'log';
    ax.XLim = [ds(1) ds(end)];
    ax.Box = 'off';
    ax = subplot(1,2,2);
    errorbar(ds,mean(dprime),std(dprime)/sqrt(ntrials),'LineWidth',2)
    xlabel('d')
    ylabel('d''')
    title(sprintf('p = %d',p))
    ax.XScale = 'log';
    ax.YScale = 'log';
    ax.XLim = [ds(1) ds(end)];
    ax.Box = 'off';
    save_fig(gcf,[rerfPath 'RandomerForest/Figures/Random_matrix_removed_columns_p_' num2str(p)])
    
    % log10(dprime) = slope*log10(d), where slope is a function of p
    slope(i) = (log10(mean(dprime(:,end))) - log10(mean(dprime(:,1))))/(log10(ds(end)) - log10(ds(1)));
end

figure;
plot(ps,slope,'LineWidth',2)
xlabel('p')
ylabel('slope')
save_fig(gcf,[rerfPath 'RandomerForest/Figures/Random_matrix_dprime_vs_d_slope'])

save([rerfPath 'RandomerForest/Results/Random_matrix_adjustment_factor'],'ps','slope')