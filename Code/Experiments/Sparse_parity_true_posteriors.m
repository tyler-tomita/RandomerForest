%% Compute Sparse Parity True Class Posteriors

close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

d = 20;
dgood = 3;
NWorkers = 2;
if d <= 5
    mtrys = 1:d;
else
    mtrys = ceil(d.^[0 1/4 1/2 3/4 1]);
end
Sigma = 1/32*ones(1,d);

xmin = -.5;
xmax = 1.5;
ymin = xmin;
ymax = xmax;
npoints = 50;
[x1gv,x2gv] = meshgrid(linspace(xmin,xmax,npoints),linspace(ymin,ymax,npoints));
X1post = x1gv(:);
X2post = x2gv(:);
X3post = zeros(npoints^2,d-2);
Xpost = [X1post X2post X3post];
n = size(Xpost,1);

allClusters = double(all_binary_sets(d));
nClust = size(allClusters,1);
inc = round(nClust*.01);
for i = 1:nClust
    if mod(i,inc)==0
        fprintf('%d %% complete\n',i/inc)
    end
    f_k(:,i) = mvnpdf(Xpost,allClusters(i,:),Sigma);
end

g = sum(f_k,2);

truth.posteriors(:,2) = sum(f_k(:,mod(sum(allClusters(:,1:dgood),2),2)==1),2)./g;
truth.posteriors(:,1) = 1 - truth.posteriors(:,2);

save([rerfPath 'RandomerForest/Results/Sparse_parity_true_posteriors.mat'],'d','dgood','Sigma','Xpost','truth')