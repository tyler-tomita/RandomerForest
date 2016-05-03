close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

d = 100;
k = 1000;

M = randmat(d,k,'sparse',1/d);
nnzs = sum(M~=0);
mn = min(nnzs);
mx = max(nnzs);

h = histogram(nnzs,'Normalization','probability');
xlabel('# nonzeros per column of A')
ylabel('Relative Frequency')
title('Sparse Rademacher Matrix')
ax = gca;
ax.Box = 'off';
ax.XTick = mn:mx;
ax.FontSize = 20;

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Rademacher_histogram'])