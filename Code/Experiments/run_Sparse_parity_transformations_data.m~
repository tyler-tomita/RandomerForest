close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

n = 1000;
dims = [2 5 10 25 50 100];
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
    dgood = min(3,d);
    x = zeros(n,d,ntrials);
    x_rot = zeros(n,d,ntrials);
    x_scale = zeros(n,d,ntrials);
    x_affine = zeros(n,d,ntrials);
    y = cell(n,ntrials);
    Sigma = 1/32*ones(1,d);
    for trial = 1:ntrials
        Mu = sparse(n,d);
        for jj = 1:n
            Mu(jj,:) = binornd(1,0.5,1,d);
            x(jj,:,trial) = mvnrnd(Mu(jj,:),Sigma);
        end
        nones = sum(Mu(:,1:dgood),2);
        y(:,trial) = cellstr(num2str(mod(nones,2)));
        R = random_rotation(d);
        S = random_scaling(n,d,0,10);
        Sigma_outlier = 4*Sigma;
        x_rot(:,:,trial) = x(:,:,trial)*R;
        x_scale(:,:,trial) = x(:,:,trial).*S;
        x_affine(:,:,trial) = (x(:,:,trial)*R).*S;
        n_out = round(0.2*n);
        xx_out = zeros(n_out,d);
        Mu_out = zeros(n_out,d);
        for jj = 1:n_out
            Mu_out(jj,:) = binornd(1,0.5,1,d);
            xx_out(jj,:) = mvnrnd(Mu_out(jj,:),Sigma_outlier);
        end
        nones_out = sum(Mu_out(:,1:dgood),2);
        yy_out = cellstr(num2str(mod(nones_out,2)));
        x_out(:,:,trial) = cat(1,x(:,:,trial),xx_out);
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

save([rerfPath 'RandomerForest/Results/Sparse_parity_transformations_data.mat'],'X','X_rot','X_scale','X_affine','X_out','Y','Y_out','n','dims','ntrials')