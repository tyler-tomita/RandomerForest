close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

n = 100;
dims = [2 10 50 100 500];
ndims = length(dims);
ntrials = 10;
Class = [0;1];
X = cell(1,ndims);
X_rot = cell(1,ndims);
X_scale = cell(1,ndims);
X_affine = cell(1,ndims);
X_out = cell(1,ndims);
Y = cell(1,ndims);
Y_out = cell(1,ndims);

for i = 1:ndims
    d = dims(i);
    d_idx = 1:d;
    mu1 = 1./sqrt(d_idx);
    mu0 = -1*mu1;
    Mu = cat(1,mu0,mu1);
    Sigma = ones(1,d);
    obj = gmdistribution(Mu,Sigma);
    x = zeros(n,d,ntrials);
    x_rot = zeros(n,d,ntrials);
    x_scale = zeros(n,d,ntrials);
    x_affine = zeros(n,d,ntrials);
    y = cell(n,ntrials);
    for trial = 1:ntrials
        [x(:,:,trial),idx] = random(obj,n);
        y(:,trial) = cellstr(num2str(Class(idx)));
        R = random_rotation(d);
        S = 10.^random_scaling(n,d,-5,5);
        Sigma_outlier = 16*Sigma;
        x_rot(:,:,trial) = x(:,:,trial)*R;
        x_scale(:,:,trial) = x(:,:,trial).*S;
        x_affine(:,:,trial) = (x(:,:,trial)*R).*S;
        outlier_model = gmdistribution(Mu,Sigma_outlier);
        [xx_out,idx_out] = random(outlier_model,0.05*n);
        x_out(:,:,trial) = cat(1,x(:,:,trial),xx_out);
        yy_out = cellstr(num2str(Class(idx_out)));
        y_out(:,trial) = cat(1,y(:,trial),yy_out);
    end
    X{i} = x;
    X_rot{i} = x_rot;
    X_scale{i} = x_scale;
    X_affine{i} = x_affine;
    X_out{i} = x_out;
    Y{i} = y;
    Y_out{i} = y_out;
    clear x_out y_out
end

save([rerfPath 'RandomerForest/Results/Trunk_transformations_data.mat'],'X','X_rot','X_scale','X_affine','X_out','Y','Y_out','n','dims','ntrials')
