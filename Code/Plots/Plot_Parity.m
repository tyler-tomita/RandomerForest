clear
close all
clc

runSims = false;

if runSims
    run_Parity
else
    load Parity.mat
end

[Lhat.rf,minIdx.rf] = min(mean_err_rf(end,:,:),[],2);
[Lhat.rerf,minIdx.rerf] = min(mean_err_rerf(end,:,:),[],2);
[Lhat.rerfdn,minIdx.rerfdn] = min(mean_err_rerfdn(end,:,:),[],2);
[Lhat.rf_rot,minIdx.rf_rot] = min(mean_err_rf_rot(end,:,:),[],2);

for i = 1:length(dims)
    sem.rf(i) = sem_rf(end,minIdx.rf(i),i);
    sem.rerf(i) = sem_rerf(end,minIdx.rerf(i),i);
    sem.rerfdn(i) = sem_rerfdn(end,minIdx.rerfdn(i),i);
    sem.rf_rot(i) = sem_rf_rot(end,minIdx.rf_rot(i),i);
end

classifiers = fieldnames(Lhat);

for i = 1:length(classifiers)
    cl = classifiers{i};
    h(i) = errorbar(dims,Lhat.(cl)(:)',sem.(cl));
    i = i + 1;
    hold on
end
xlabel('d')
ylabel('Lhat')
title('Parity')
legend('RF','RerF','RerFdn','Rotation RF')

save_fig(gcf,'~/RandomerForest/Figures/Parity_2015_11_9')