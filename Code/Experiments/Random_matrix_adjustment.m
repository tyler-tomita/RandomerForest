%% Plot number of unique columns vs number of total columns using adjustment factors

close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(123);

load Random_matrix_adjustment_factor

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
    
    dprime = zeros(ntrials,length(ds));
    
    for j = 1:length(ds)
        d = ds(j);
        for trial = 1:ntrials
            fprintf('trial %d\n',trial)
            M = randmat(p,d,'sparse-adjusted',1/p,ceil(d^(1/interp1(ps,slope,p))));
            dprime(trial,j) = size(M,2);
        end
    end
    figure(i)
    ax = gca;
    errorbar(ds,mean(dprime),std(dprime)/sqrt(ntrials),'LineWidth',2)
    xlabel('d')
    ylabel('adjusted d''')
    title(sprintf('p = %d',p))
%     ax.XScale = 'log';
%     ax.YScale = 'log';
    ax.XLim = [ds(1) ds(end)];
    ax.Box = 'off';
    save_fig(gcf,[rerfPath 'RandomerForest/Figures/Random_matrix_adjustment_p_' num2str(p)])
end