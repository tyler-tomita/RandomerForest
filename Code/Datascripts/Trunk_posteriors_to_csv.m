clear
close all
clc

load Trunk_vary_n_data

for j = 1:length(ps)
    p = ps(j);
    dlmwrite(sprintf('~/Documents/R/Data/Trunk/dat/Test/Trunk_test_set_posteriors_p%d.dat',p),...
    ClassPosteriors{j},...
    'delimiter','\t','precision','%0.15f');
end