clear
close all
clc

load Sparse_parity_data

for i = 1:length(Xtrain)
    p = ps(i);
    for trial = 1:ntrials
        dlmwrite(sprintf('~/Documents/R/Data/Sparse_parity/dat/Train/Sparse_parity_train_set_p%d_trial%d.dat',p,trial),...
            [Xtrain(i).Untransformed(:,:,trial),cellfun(@str2num,Ytrain(i).Untransformed(:,trial))],...
            'delimiter','\t','precision','%0.15f');
    end
    dlmwrite(sprintf('~/Documents/R/Data/Sparse_parity/dat/Test/Sparse_parity_test_set_p%d.dat',p),...
    [Xtest(i).Untransformed(:,:,1),cellfun(@str2num,Ytest(i).Untransformed(:,1))],...
    'delimiter','\t','precision','%0.15f');
end