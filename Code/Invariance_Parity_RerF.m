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
%dims = round(logspace(log10(2),3,10));
dims = [2 3 4 5 10];
ntrees = 1000;
ntrials = 10;
NWorkers = 2;
Class = [0;1];

rerf_err = NaN(ntrials,length(dims));
rerf_rot_err = NaN(ntrials,length(dims));
rerf_trans_err = NaN(ntrials,length(dims));
rerf_scale_err = NaN(ntrials,length(dims));
rerf_affine_err = NaN(ntrials,length(dims));
rerf_out_err = NaN(ntrials,length(dims));

for trial = 1:ntrials
    
    fprintf('trial %d\n',trial)
    
    for i = 1:length(dims)
    
        d = dims(i);
        fprintf('d = %d\n',d)
        nvartosample = ceil(d^(2/3));
        d_idx = 1:d;
        X = zeros(n,d);
        Sigma = 1/32*ones(1,d);
        Mu = zeros(n,d);
        for j = 1:n
            Mu(j,:) = binornd(1,0.5,1,d);
            X(j,:) = mvnrnd(Mu(j,:),Sigma);
        end
        
        nones = sum(Mu,2);
        Y = cellstr(num2str(mod(nones,2)));
        
        R = random_rotation(d);
        T = random_translation(n,d,-1,1);
        S = random_scaling(n,d,0,10);
        Sigma_outlier = 4*Sigma;
        X_rot = X*R;
        X_trans = X + T;
        X_scale = X.*S;
        X_affine = (X*R).*S;
        outlier_model = gmdistribution(Mu,Sigma_outlier);
        n_out = round(0.2*n);
        X_out = zeros(n_out,d);
        Mu_out = zeros(n_out,d);
        for j = 1:n_out
            Mu_out(j,:) = binornd(1,0.5,1,d);
            X_out(j,:) = mvnrnd(Mu(j,:),Sigma);
        end
        nones_out = sum(Mu_out,2);
        Y_out = cellstr(num2str(mod(nones_out,2)));
        X_out = cat(1,X,X_out);
        Y_out = cat(1,Y,Y_out);

        rerf = rpclassificationforest(ntrees,X,Y,'nvartosample',nvartosample,'mdiff','off','sparsemethod','jovo','Robust',false,'NWorkers',NWorkers);
        rerf_err(trial,i) = oobpredict(rerf,X,Y,'last');
        
        rerf_rot = rpclassificationforest(ntrees,X_rot,Y,'nvartosample',nvartosample,'mdiff','off','sparsemethod','jovo','Robust',false,'NWorkers',NWorkers);
        rerf_rot_err(trial,i) = oobpredict(rerf_rot,X_rot,Y,'last');
        clear rerf_rot

        rerf_trans = rpclassificationforest(ntrees,X_trans,Y,'nvartosample',nvartosample,'mdiff','off','sparsemethod','jovo','Robust',false,'NWorkers',NWorkers);
        rerf_trans_err(trial,i) = oobpredict(rerf_trans,X_trans,Y,'last');
        clear rerf_trans

        rerf_scale = rpclassificationforest(ntrees,X_scale,Y,'nvartosample',nvartosample,'mdiff','off','sparsemethod','jovo','Robust',false,'NWorkers',NWorkers);
        rerf_scale_err(trial,i) = oobpredict(rerf_scale,X_scale,Y,'last');
        clear rerf_scale
        
        rerf_affine = rpclassificationforest(ntrees,X_affine,Y,'nvartosample',nvartosample,'mdiff','off','sparsemethod','jovo','Robust',false,'NWorkers',NWorkers);
        rerf_affine_err(trial,i) = oobpredict(rerf_affine,X_affine,Y,'last');
        clear rerf_affine

        rerf_out = rpclassificationforest(ntrees,X_out,Y_out,'nvartosample',nvartosample,'mdiff','off','sparsemethod','jovo','Robust',false,'NWorkers',NWorkers);
        rerf_out_err(trial,i) = oobpredict(rerf_out,X_out,Y_out,'last');
        clear rerf_out
    end
end

mean_rerf_err = mean(rerf_err);
mean_rerf_rot_err = mean(rerf_rot_err);
mean_rerf_trans_err = mean(rerf_trans_err);
mean_rerf_scale_err = mean(rerf_scale_err);
mean_rerf_affine_err = mean(rerf_affine_err);
mean_rerf_out_err = mean(rerf_out_err);

sem_rerf = std(rerf_err)/sqrt(ntrials);
sem_rerf_rot = std(rerf_rot_err)/sqrt(ntrials);
sem_rerf_trans = std(rerf_trans_err)/sqrt(ntrials);
sem_rerf_scale = std(rerf_scale_err)/sqrt(ntrials);
sem_rerf_affine = std(rerf_affine_err)/sqrt(ntrials);
sem_rerf_out = std(rerf_out_err)/sqrt(ntrials);

Ynames = {'mean_rerf_err' 'mean_rerf_rot_err' 'mean_rerf_trans_err' 'mean_rerf_scale_err' 'mean_rerf_affine_err' 'mean_rerf_out_err'};
Enames = {'sem_rerf' 'sem_rerf_rot' 'sem_rerf_trans' 'sem_rerf_scale' 'sem_rerf_affine' 'sem_rerf_out'};
lspec = {'-bs','-rs','-gs' '-cs' '-ms' '-ks'};
facespec = {'b','r','g' 'c' 'm' 'k'};
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
    hold on
end
set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel(sprintf('OOB Error for %d Trees',ntrees))
legend('Untransformed','Rotated','Translated','Scaled','Affine','Outlier')
title('RerF')

filename = 'Invariance_Parity_RerF';
save_fig(gcf,filename)
