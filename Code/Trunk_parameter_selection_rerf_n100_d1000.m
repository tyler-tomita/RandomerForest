close all
clear
clc

n = 100;
d =1000;
ntrees = 1000;
ntrials = 20;
NWorkers = 2;
Colors = linspecer(8,'sequential');
mtrys = ceil(d.^[1/4 1/2 3/4 1]);
err_rf = zeros(ntrees,length(mtrys),ntrials);
err_rerf = zeros(ntrees,length(mtrys),ntrials);
Class = [0;1];

for trial = 1:ntrials
    
    fprintf('trial %d\n',trial)
    
    d_idx = 1:d;
    mu1 = 1./sqrt(d_idx);
    mu0 = -1*mu1;
    Mu = cat(1,mu0,mu1);
    Sigma = ones(1,d);
    obj = gmdistribution(Mu,Sigma);
    [X,idx] = random(obj,n);
    Ystr = cellstr(num2str(Class(idx)));
    
    i = 1;
    
    for mtry = mtrys
        
        fprintf('mtry = %d\n',mtry)

        rf = rpclassificationforest(ntrees,X,Ystr,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers);
        err_rf(:,i,trial) = oobpredict(rf,X,Ystr,'every');
        
        rerf = rpclassificationforest(ntrees,X,Ystr,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers);
        err_rerf(:,i,trial) = oobpredict(rerf,X,Ystr,'every');
        
        i = i + 1;
    end
end

save(sprintf('Trunk_parameter_selection_rerf_n%d_d%d.mat',n,d),'err_rf','err_rerf')

mean_err_rf = mean(err_rf,3);
mean_err_rerf = mean(err_rerf,3);

i = 1;

for mtry = mtrys
    plot(mean_err_rf(:,i),'Color',Colors(i,:))
    hold on
    i = i + 1;
end

for mtry = mtrys
    plot(mean_err_rerf(:,i-4),'Color',Colors(i,:))
    hold on
    i = i + 1;
end

legend('rf (mtry = d^1^/^4)','rf (mtry = d^1^/^2)','rf (mtry = d^3^/^4)','rf (mtry = d)','rerf (mtry = d^1^/^4)','rerf (mtry = d^1^/^2)','rerf (mtry = d^3^/^4)','rerf (mtry = d)')

xlabel('# trees')
ylabel('oob error')
title(sprintf('Trunk (n=%d, d=%d, ntrials=%d)',n,d,ntrials))

save_fig(gcf,sprintf('Trunk_parameter_selection_rerf_n%d_d%d',n,d))
