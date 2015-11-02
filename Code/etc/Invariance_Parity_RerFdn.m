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

rerfdn_err = NaN(ntrials,length(dims));
rerfdn_rot_err = NaN(ntrials,length(dims));
rerfdn_trans_err = NaN(ntrials,length(dims));
rerfdn_scale_err = NaN(ntrials,length(dims));
rerfdn_affine_err = NaN(ntrials,length(dims));
rerfdn_out_err = NaN(ntrials,length(dims));

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

        rerfdn = rpclassificationforest(ntrees,X,Y,'nvartosample',nvartosample,'mdiff','node','sparsemethod','sparse','Robust',false,'NWorkers',NWorkers);
        rerfdn_err(trial,i) = oobpredict(rerfdn,X,Y,'last');
        
        rerfdn_rot = rpclassificationforest(ntrees,X_rot,Y,'nvartosample',nvartosample,'mdiff','node','sparsemethod','sparse','Robust',false,'NWorkers',NWorkers);
        rerfdn_rot_err(trial,i) = oobpredict(rerfdn_rot,X_rot,Y,'last');
        clear rerfdn_rot

        rerfdn_trans = rpclassificationforest(ntrees,X_trans,Y,'nvartosample',nvartosample,'mdiff','node','sparsemethod','sparse','Robust',false,'NWorkers',NWorkers);
        rerfdn_trans_err(trial,i) = oobpredict(rerfdn_trans,X_trans,Y,'last');
        clear rerfdn_trans

        rerfdn_scale = rpclassificationforest(ntrees,X_scale,Y,'nvartosample',nvartosample,'mdiff','node','sparsemethod','sparse','Robust',false,'NWorkers',NWorkers);
        rerfdn_scale_err(trial,i) = oobpredict(rerfdn_scale,X_scale,Y,'last');
        clear rerfdn_scale
        
        rerfdn_affine = rpclassificationforest(ntrees,X_affine,Y,'nvartosample',nvartosample,'mdiff','node','sparsemethod','sparse','Robust',false,'NWorkers',NWorkers);
        rerfdn_affine_err(trial,i) = oobpredict(rerfdn_affine,X_affine,Y,'last');
        clear rerfdn_affine

        rerfdn_out = rpclassificationforest(ntrees,X_out,Y_out,'nvartosample',nvartosample,'mdiff','node','sparsemethod','sparse','Robust',false,'NWorkers',NWorkers);
        rerfdn_out_err(trial,i) = oobpredict(rerfdn_out,X_out,Y_out,'last');
        clear rerfdn_out
    end
end

mean_rerfdn_err = mean(rerfdn_err);
mean_rerfdn_rot_err = mean(rerfdn_rot_err);
mean_rerfdn_trans_err = mean(rerfdn_trans_err);
mean_rerfdn_scale_err = mean(rerfdn_scale_err);
mean_rerfdn_affine_err = mean(rerfdn_affine_err);
mean_rerfdn_out_err = mean(rerfdn_out_err);

sem_rerfdn = std(rerfdn_err)/sqrt(ntrials);
sem_rerfdn_rot = std(rerfdn_rot_err)/sqrt(ntrials);
sem_rerfdn_trans = std(rerfdn_trans_err)/sqrt(ntrials);
sem_rerfdn_scale = std(rerfdn_scale_err)/sqrt(ntrials);
sem_rerfdn_affine = std(rerfdn_affine_err)/sqrt(ntrials);
sem_rerfdn_out = std(rerfdn_out_err)/sqrt(ntrials);

Ynames = {'mean_rerfdn_err' 'mean_rerfdn_rot_err' 'mean_rerfdn_trans_err' 'mean_rerfdn_scale_err' 'mean_rerfdn_affine_err' 'mean_rerfdn_out_err'};
Enames = {'sem_rerfdn' 'sem_rerfdn_rot' 'sem_rerfdn_trans' 'sem_rerfdn_scale' 'sem_rerfdn_affine' 'sem_rerfdn_out'};
lspec = {'-bs','-rs','-gs' '-cs' '-ms' '-ks'};
facespec = {'b','r','g' 'c' 'm' 'k'};
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
    hold on
end
xlabel('# Ambient Dimensions')
ylabel(sprintf('OOB Error for %d Trees',ntrees))
legend('Untransformed','Rotated','Translated','Scaled','Affine','Outlier')
title('rerfdn')

filename = 'Invariance_Parity_rerfdn';
save_fig(gcf,filename)
