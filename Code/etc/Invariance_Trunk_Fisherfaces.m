close all
clear
clc

fpath = mfilename('fullpath');
findex = strfind(fpath,'/');
rootDir=fpath(1:findex(end-1));
p = genpath(rootDir);
gits=strfind(p,'.git');
colons=strfind(p,':');
for i=0:length(gits)-1
endGit=find(colons>gits(end-i),1);
p(colons(endGit-1):colons(endGit)-1)=[];
end
addpath(p);

n = 100;    %# samples
dims = round(logspace(log10(2),3,10));
ntrials = 10;
Class = [0;1];

ff_err = NaN(ntrials,length(dims));
ff_rot_err = NaN(ntrials,length(dims));
ff_trans_err = NaN(ntrials,length(dims));
ff_scale_err = NaN(ntrials,length(dims));
ff_affine_err = NaN(ntrials,length(dims));
ff_out_err = NaN(ntrials,length(dims));

for trial = 1:ntrials
    
    fprintf('trial %d\n',trial)
    
    for i = 1:length(dims)
    
        d = dims(i);
        fprintf('d = %d\n',d)
        nvartosample = ceil(d^(2/3));
        d_idx = 1:d;
        mu1 = 1./sqrt(d_idx);
        mu0 = -1*mu1;
        Mu = cat(1,mu0,mu1);
        Sigma = 1*speye(d);
        obj = gmdistribution(Mu,Sigma);
        [X,idx] = random(obj,n);
        Y = cellstr(num2str(Class(idx)));
        
        R = random_rotation(d);
        T = random_translation(n,d,-1,1);
        S = random_scaling(n,d,0,10);
        Sigma_outlier = 16*Sigma;
        X_rot = X*R;
        X_trans = X + T;
        X_scale = X.*S;
        X_affine = (X*R).*S;
        outlier_model = gmdistribution(Mu,Sigma_outlier);
        [X_out,idx_out] = random(outlier_model,0.05*n);
        X_out = cat(1,X,X_out);
        Y_out = cellstr(num2str(Class(idx_out)));
        Y_out = cat(1,Y,Y_out);

        ff_err(trial,i) = fisherface_loocv(X,Y);
        ff_rot_err(trial,i) = fisherface_loocv(X_rot,Y);
        ff_trans_err(trial,i) = fisherface_loocv(X_trans,Y);
        ff_scale_err(trial,i) = fisherface_loocv(X_scale,Y);
        ff_affine_err(trial,i) = fisherface_loocv(X_affine,Y);
        ff_out_err(trial,i) = fisherface_loocv(X_out,Y_out);
    end
end

save('Invariance_Trunk_Fisherfaces.mat','ff_err','ff_rot_err',...
    'ff_trans_err','ff_scale_err','ff_affine_err','ff_out_err')

mean_ff_err = mean(ff_err);
mean_ff_rot_err = mean(ff_rot_err);
mean_ff_trans_err = mean(ff_trans_err);
mean_ff_scale_err = mean(ff_scale_err);
mean_ff_affine_err = mean(ff_affine_err);
mean_ff_out_err = mean(ff_out_err);

sem_ff = std(ff_err)/sqrt(ntrials);
sem_ff_rot = std(ff_rot_err)/sqrt(ntrials);
sem_ff_trans = std(ff_trans_err)/sqrt(ntrials);
sem_ff_scale = std(ff_scale_err)/sqrt(ntrials);
sem_ff_affine = std(ff_affine_err)/sqrt(ntrials);
sem_ff_out = std(ff_out_err)/sqrt(ntrials);

Ynames = {'mean_ff_err' 'mean_ff_rot_err' 'mean_ff_trans_err' 'mean_ff_scale_err' 'mean_ff_affine_err' 'mean_ff_out_err'};
Enames = {'sem_ff' 'sem_ff_rot' 'sem_ff_trans' 'sem_ff_scale' 'sem_ff_affine' 'sem_ff_out'};
lspec = {'-bs','-rs','-gs' '-cs' '-ms' '-ks'};
facespec = {'b','r','g' 'c' 'm' 'k'};
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
    hold on
end
set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel('Lhat')
legend('Untransformed','Rotated','Translated','Scaled','Affine','Outlier')
title('Fisherfaces')


filename = 'Invariance_Trunk_Fisherfaces';
save_fig(gcf,filename)