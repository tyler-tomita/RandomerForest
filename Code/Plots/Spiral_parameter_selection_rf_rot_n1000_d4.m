close all
clear
clc

n = 1000;
d = 4;
[X,Y] = Spiral_2class(n,d);
Ystr = cellstr(num2str(Y));
ntrees = 500;
ntrials = 20;
NWorkers = 2;
mtrys = 1:d;
err_rf = zeros(ntrees,length(mtrys),ntrials);
err_rf_rot = zeros(ntrees,length(mtrys),ntrials);
Colors = linspecer(2*length(mtrys),'sequential');

for trial = 1:ntrials
    
    fprintf('trial %d\n',trial)
    
    i = 1;
    
    for mtry = mtrys
        
        fprintf('mtry = %d\n',mtry)

        rf = rpclassificationforest(ntrees,X,Ystr,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers);
        err_rf(:,i,trial) = oobpredict(rf,X,Ystr,'every');
        
        rf_rot = rpclassificationforest(ntrees,X,Ystr,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'rotate',true);
        err_rf_rot(:,i,trial) = oobpredict(rf_rot,X,Ystr,'every');
        
        i = i + 1;
    end
end

save(sprintf('Spiral_parameter_selection_rf_rot_n%d_d%d.mat',n,d),'err_rf','err_rf_rot')

sem_rf = std(err_rf,[],3)/sqrt(ntrials);
sem_rf_rot = std(err_rf_rot,[],3)/sqrt(ntrials);

var_rf = var(err_rf,0,3);
var_rf_rot = var(err_rf_rot,0,3);

mean_err_rf = mean(err_rf,3);
mean_err_rf_rot = mean(err_rf_rot,3);

legend_names = {};

i = 1;

for mtry = mtrys
    h(i) = plot(mean_err_rf(:,i),'Color',Colors(i,:));
    hold on
    legend_names{i} = sprintf('rf (mtry = %d)',mtry);
    i = i + 1;
end

for mtry = mtrys
    h(i) = plot(mean_err_rf_rot(:,i-length(mtrys)),'Color',Colors(i,:));
    hold on
    legend_names{i} = sprintf('rf rotate (mtry = %d)',mtry);
    i = i + 1;
end

legend(h,legend_names)
xlabel('# trees')
ylabel('oob error')
title(sprintf('Spiral (n=%d, d=%d, ntrials=%d)',n,d,ntrials))

save_fig(gcf,sprintf('Spiral_parameter_selection_rf_rot_n%d_d%d',n,d))

figure(2)

i = 1;

for mtry = mtrys
    h(i) = plot(var_rf(:,i),'Color',Colors(i,:));
    hold on
    i = i + 1;
end

for mtry = mtrys
    h(i) = plot(var_rf_rot(:,i-length(mtrys)),'Color',Colors(i,:));
    hold on
    i = i + 1;
end

legend(h,legend_names)
xlabel('# trees')
ylabel('var(oob error)')
title(sprintf('Spiral (n=%d, d=%d, ntrials=%d)',n,d,ntrials))

save_fig(gcf,sprintf('Spiral_parameter_selection_rf_rot_variance_n%d_d%d',n,d))

figure(3)

i = 1;

for mtry = mtrys
    h(i) = errorbar(1:ntrees,mean_err_rf(:,i),sem_rf(:,i),'Color',Colors(i,:));
    hold on
    i = i + 1;
end

for mtry = mtrys
    h(i) = errorbar(1:ntrees,mean_err_rf_rot(:,i-length(mtrys)),sem_rf_rot(:,i-length(mtrys)),'Color',Colors(i,:));
    hold on
    i = i + 1;
end

legend(h,legend_names)
xlabel('# trees')
ylabel('oob error')
title(sprintf('Spiral (n=%d, d=%d, ntrials=%d)',n,d,ntrials))

save_fig(gcf,sprintf('Spiral_parameter_selection_rf_rotu_errorbar_n%d_d%d',n,d))