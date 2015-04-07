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
dims = round(logspace(log10(2),3,7));
ntrees = 1000;
ntrials = 10;
Class = [0;1];

rf_err = NaN(ntrials,length(dims));
f1_err = NaN(ntrials,length(dims));
f2_err = NaN(ntrials,length(dims));
rf_rot_err = NaN(ntrials,length(dims));
f1_rot_err = NaN(ntrials,length(dims));
f2_rot_err = NaN(ntrials,length(dims));
rf_trans_err = NaN(ntrials,length(dims));
f1_trans_err = NaN(ntrials,length(dims));
f2_trans_err = NaN(ntrials,length(dims));
rf_scale_err = NaN(ntrials,length(dims));
f1_scale_err = NaN(ntrials,length(dims));
f2_scale_err = NaN(ntrials,length(dims));
rf_affine_err = NaN(ntrials,length(dims));
f1_affine_err = NaN(ntrials,length(dims));
f2_affine_err = NaN(ntrials,length(dims));
rf_out_err = NaN(ntrials,length(dims));
f1_out_err = NaN(ntrials,length(dims));
f2_out_err = NaN(ntrials,length(dims));

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
        Sigma_outlier = 4*Sigma;
        X_rot = X*R;
        X_trans = X + T;
        X_scale = X.*S;
        X_affine = (X*R).*S;
        outlier_model = gmdistribution(Mu,Sigma_outlier);
        [X_out,idx_out] = random(outlier_model,0.2*n);
        X_out = cat(1,X,X_out);
        Y_out = cellstr(num2str(Class(idx_out)));
        Y_out = cat(1,Y,Y_out);

        rf = rpclassificationforest(ntrees,X,Y,'nvartosample',nvartosample,'mdiff','off','RandomForest',true);
        rf_err(trial,i) = oobpredict(rf,X,Y,'last');
        clear rf

        f1 = rpclassificationforest(ntrees,X,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        f1_err(trial,i) = oobpredict(f1,X,Y,'last');
        clear f1
        
        f2 = rpclassificationforest(ntrees,X,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new','Robust',true);
        f2_err(trial,i) = oobpredict(f2,X,Y,'last');
        clear f2

        rf_rot = rpclassificationforest(ntrees,X_rot,Y,'nvartosample',nvartosample,'mdiff','off','RandomForest',true);
        rf_rot_err(trial,i) = oobpredict(rf_rot,X_rot,Y,'last');
        clear rf_rot

        f1_rot = rpclassificationforest(ntrees,X_rot,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        f1_rot_err(trial,i) = oobpredict(f1_rot,X_rot,Y,'last');
        clear f1_rot
        
        f2_rot = rpclassificationforest(ntrees,X_rot,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new','Robust',true);
        f2_rot_err(trial,i) = oobpredict(f2_rot,X_rot,Y,'last');
        clear f2_rot

        rf_trans = rpclassificationforest(ntrees,X_trans,Y,'nvartosample',nvartosample,'mdiff','off','RandomForest',true);
        rf_trans_err(trial,i) = oobpredict(rf_trans,X_trans,Y,'last');
        clear rf_trans

        f1_trans = rpclassificationforest(ntrees,X_trans,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        f1_trans_err(trial,i) = oobpredict(f1_trans,X_trans,Y,'last');
        clear f1_trans
        
        f2_trans = rpclassificationforest(ntrees,X_trans,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new','Robust',true);
        f2_trans_err(trial,i) = oobpredict(f2_trans,X_trans,Y,'last');
        clear f2_trans

        rf_scale = rpclassificationforest(ntrees,X_scale,Y,'nvartosample',nvartosample,'mdiff','off','RandomForest',true);
        rf_scale_err(trial,i) = oobpredict(rf_scale,X_scale,Y,'last');
        clear rf_scale

        f1_scale = rpclassificationforest(ntrees,X_scale,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        f1_scale_err(trial,i) = oobpredict(f1_scale,X_scale,Y,'last');
        clear f1_scale
        
        f2_scale = rpclassificationforest(ntrees,X_scale,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new','Robust',true);
        f2_scale_err(trial,i) = oobpredict(f2_scale,X_scale,Y,'last');
        clear f2_scale
        
        rf_affine = rpclassificationforest(ntrees,X_affine,Y,'nvartosample',nvartosample,'mdiff','off','RandomForest',true);
        rf_affine_err(trial,i) = oobpredict(rf_affine,X_affine,Y,'last');
        clear rf_affine

        f1_affine = rpclassificationforest(ntrees,X_affine,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        f1_affine_err(trial,i) = oobpredict(f1_affine,X_affine,Y,'last');
        clear f1_affine
        
        f2_affine = rpclassificationforest(ntrees,X_affine,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new','Robust',true);
        f2_affine_err(trial,i) = oobpredict(f2_affine,X_affine,Y,'last');
        clear f2_affine

        
        rf_out = rpclassificationforest(ntrees,X_out,Y_out,'nvartosample',nvartosample,'mdiff','off','RandomForest',true);
        rf_out_err(trial,i) = oobpredict(rf_out,X_out,Y_out,'last');
        clear rf_out

        f1_out = rpclassificationforest(ntrees,X_out,Y_out,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        f1_out_err(trial,i) = oobpredict(f1_out,X_out,Y_out,'last');
        clear f1_out
        
        f2_out = rpclassificationforest(ntrees,X_out,Y_out,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new','Robust',true);
        f2_out_err(trial,i) = oobpredict(f2_out,X_out,Y_out,'last');
        clear f2_out
    end
end

save('Invariance_Trunk.mat','R','rf_err','f1_err','f2_err','rf_rot_err',...
    'f1_rot_err','f2_rot_err','rf_trans_err','f1_trans_err','f2_trans_err',...
    'rf_scale_err','f1_scale_err','f2_scale_err','rf_affine_err','f1_affine_err','f2_affine_err',...
    'rf_out_err','f1_out_err','f2_out_err')

mean_rf_err = mean(rf_err);
mean_f1_err = mean(f1_err);
mean_f2_err = mean(f2_err);
mean_rf_rot_err = mean(rf_rot_err);
mean_f1_rot_err = mean(f1_rot_err);
mean_f2_rot_err = mean(f2_rot_err);
mean_rf_trans_err = mean(rf_trans_err);
mean_f1_trans_err = mean(f1_trans_err);
mean_f2_trans_err = mean(f2_trans_err);
mean_rf_scale_err = mean(rf_scale_err);
mean_f1_scale_err = mean(f1_scale_err);
mean_f2_scale_err = mean(f2_scale_err);
mean_rf_affine_err = mean(rf_affine_err);
mean_f1_affine_err = mean(f1_affine_err);
mean_f2_affine_err = mean(f2_affine_err);
mean_rf_out_err = mean(rf_out_err);
mean_f1_out_err = mean(f1_out_err);
mean_f2_out_err = mean(f2_out_err);

sem_rf = std(rf_err)/sqrt(ntrials);
sem_f1 = std(f1_err)/sqrt(ntrials);
sem_f2 = std(f2_err)/sqrt(ntrials);
sem_rf_rot = std(rf_rot_err)/sqrt(ntrials);
sem_f1_rot = std(f1_rot_err)/sqrt(ntrials);
sem_f2_rot = std(f2_rot_err)/sqrt(ntrials);
sem_rf_trans = std(rf_trans_err)/sqrt(ntrials);
sem_f1_trans = std(f1_trans_err)/sqrt(ntrials);
sem_f2_trans = std(f2_trans_err)/sqrt(ntrials);
sem_rf_scale = std(rf_scale_err)/sqrt(ntrials);
sem_f1_scale = std(f1_scale_err)/sqrt(ntrials);
sem_f2_scale = std(f2_scale_err)/sqrt(ntrials);
sem_rf_affine = std(rf_affine_err)/sqrt(ntrials);
sem_f1_affine = std(f1_affine_err)/sqrt(ntrials);
sem_f2_affine = std(f2_affine_err)/sqrt(ntrials);
sem_rf_out = std(rf_out_err)/sqrt(ntrials);
sem_f1_out = std(f1_out_err)/sqrt(ntrials);
sem_f2_out = std(f2_out_err)/sqrt(ntrials);

Ynames = {'mean_rf_err' 'mean_rf_rot_err' 'mean_rf_trans_err' 'mean_rf_scale_err' 'mean_rf_affine_err' 'mean_rf_out_err'};
Enames = {'sem_rf' 'sem_rf_rot' 'sem_rf_trans' 'sem_rf_scale' 'sem_rf_affine' 'sem_rf_out'};
lspec = {'-bs','-rs','-gs' '-cs' '-ms' '-ks'};
facespec = {'b','r','g' 'c' 'm' 'k'};
subplot(1,3,1)
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
    hold on
end
set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel(sprintf('OOB Error for %d Trees',ntrees))
legend('Untransformed','Rotated','Translated','Scaled','Affine','Outlier')
title('Random Forest')

Ynames = {'mean_f1_err' 'mean_f1_rot_err' 'mean_f1_trans_err' 'mean_f1_scale_err' 'mean_f1_affine_err' 'mean_f1_out_err'};
Enames = {'sem_f1' 'sem_f1_rot' 'sem_f1_trans' 'sem_f1_scale' 'sem_f1_affine' 'sem_f1_out'};
subplot(1,3,2)
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
    hold on
end
set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel(sprintf('OOB Error for %d Trees',ntrees))
legend('Untransformed','Rotated','Translated','Scaled','Affine','Outlier')
title('Sparse Randomer Forest w/ Mean Diff')

Ynames = {'mean_f2_err' 'mean_f2_rot_err' 'mean_f2_trans_err' 'mean_f2_scale_err' 'mean_f2_affine_err' 'mean_f2_out_err'};
Enames = {'sem_f2' 'sem_f2_rot' 'sem_f2_trans' 'sem_f2_scale' 'sem_f2_affine' 'sem_f2_out'};
subplot(1,3,3)
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
    hold on
end
set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel(sprintf('OOB Error for %d Trees',ntrees))
legend('Untransformed','Rotated','Translated','Scaled','Affine','Outlier')
title('Robust Sparse Randomer Forest w/ Mean Diff')

filename = 'Invariance_Trunk_v2';
save_fig(gcf,filename)