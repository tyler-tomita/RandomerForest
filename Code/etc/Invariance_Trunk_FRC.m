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
nmix = 2;

frc_err = NaN(ntrials,length(dims));
frc_rot_err = NaN(ntrials,length(dims));
frc_trans_err = NaN(ntrials,length(dims));
frc_scale_err = NaN(ntrials,length(dims));
frc_affine_err = NaN(ntrials,length(dims));
frc_out_err = NaN(ntrials,length(dims));

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
        
        frc = rpclassificationforest(ntrees,X,Y,'nvartosample',nvartosample,'mdiff','off','sparsemethod','frc','Robust',false,'NWorkers',NWorkers,'nmix',nmix);
        frc_err(trial,i) = oobpredict(frc,X,Y,'last');
        
        frc_rot = rpclassificationforest(ntrees,X_rot,Y,'nvartosample',nvartosample,'mdiff','off','sparsemethod','frc','Robust',false,'NWorkers',NWorkers,'nmix',nmix);
        frc_rot_err(trial,i) = oobpredict(frc_rot,X_rot,Y,'last');
        clear frc_rot

        frc_trans = rpclassificationforest(ntrees,X_trans,Y,'nvartosample',nvartosample,'mdiff','off','sparsemethod','frc','Robust',false,'NWorkers',NWorkers,'nmix',nmix);
        frc_trans_err(trial,i) = oobpredict(frc_trans,X_trans,Y,'last');
        clear frc_trans

        frc_scale = rpclassificationforest(ntrees,X_scale,Y,'nvartosample',nvartosample,'mdiff','off','sparsemethod','frc','Robust',false,'NWorkers',NWorkers,'nmix',nmix);
        frc_scale_err(trial,i) = oobpredict(frc_scale,X_scale,Y,'last');
        clear frc_scale
        
        frc_affine = rpclassificationforest(ntrees,X_affine,Y,'nvartosample',nvartosample,'mdiff','off','sparsemethod','frc','Robust',false,'NWorkers',NWorkers,'nmix',nmix);
        frc_affine_err(trial,i) = oobpredict(frc_affine,X_affine,Y,'last');
        clear frc_affine

        frc_out = rpclassificationforest(ntrees,X_out,Y_out,'nvartosample',nvartosample,'mdiff','off','sparsemethod','frc','Robust',false,'NWorkers',NWorkers,'nmix',nmix);
        frc_out_err(trial,i) = oobpredict(frc_out,X_out,Y_out,'last');
        clear frc_out
    end
end

mean_frc_err = nanmean(frc_err);
mean_frc_rot_err = nanmean(frc_rot_err);
mean_frc_trans_err = nanmean(frc_trans_err);
mean_frc_scale_err = nanmean(frc_scale_err);
mean_frc_affine_err = nanmean(frc_affine_err);
mean_frc_out_err = nanmean(frc_out_err);

sem_frc = nanstd(frc_err)/sqrt(ntrials);
sem_frc_rot = nanstd(frc_rot_err)/sqrt(ntrials);
sem_frc_trans = nanstd(frc_trans_err)/sqrt(ntrials);
sem_frc_scale = nanstd(frc_scale_err)/sqrt(ntrials);
sem_frc_affine = nanstd(frc_affine_err)/sqrt(ntrials);
sem_frc_out = nanstd(frc_out_err)/sqrt(ntrials);

Ynames = {'mean_frc_err' 'mean_frc_rot_err' 'mean_frc_trans_err' 'mean_frc_scale_err' 'mean_frc_affine_err' 'mean_frc_out_err'};
Enames = {'sem_frc' 'sem_frc_rot' 'sem_frc_trans' 'sem_frc_scale' 'sem_frc_affine' 'sem_frc_out'};
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
title('FRC')

filename = 'Invariance_Parity_FRC';
save_fig(gcf,filename)
