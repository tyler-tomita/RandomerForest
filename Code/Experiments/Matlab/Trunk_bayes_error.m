%Analytically computes bayes error for Trunk as a function of d

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

dims = [2 10 50 100 500];
bayes_error = zeros(size(dims));

for i = 1:length(dims)
    d = dims(i);
    d_idx = 1:d;
    mu1 = 1./sqrt(d_idx);
    mu0 = -1*mu1;
    Sigma = speye(d);
    delta = mu1 - mu0;
    bayes_error(i) = 1 - normcdf(sqrt(delta*Sigma*delta')/2);
end

save([rerfPath 'RandomerForest/Results/Trunk_bayes_error.mat'],'dims','bayes_error')
