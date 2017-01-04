clear
close all
clc

load Sparse_parity_vary_n_data

for j = 1:length(ps)
    p = ps(j);
    for i = 1:length(ns{j})
        n = ns{j}(i);
        for trial = 1:ntrials
            dlmwrite(sprintf('~/Documents/R/Data/Sparse_parity/dat/Train/Sparse_parity_train_set_n%d_p%d_trial%d.dat',n,p,trial),...
                [Xtrain{i,j}(:,:,trial),cellfun(@str2num,Ytrain{i,j}(:,trial))],...
                'delimiter','\t','precision','%0.15f');
        end
    end
    dlmwrite(sprintf('~/Documents/R/Data/Sparse_parity/dat/Test/Sparse_parity_test_set_p%d.dat',p),...
    [Xtest{j},cellfun(@str2num,Ytest{j})],...
    'delimiter','\t','precision','%0.15f');
end