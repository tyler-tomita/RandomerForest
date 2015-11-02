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
dims = round(logspace(log10(2),3,5));
ntrees = 1000;
ntrials = 10;
nclasses = 4;
Classes = 1:nclasses;
ncomponents = 2;    %# of mixture components per class
J = nclasses*ncomponents; %total number of mixture components
p = ones(1,J)/J;    %mixture probabilities

%randomly allocate components evenly across all classes
Class = repmat(transpose(1:nclasses),1,ncomponents);
Class = Class(:);
Class = Class(randperm(length(Class)));

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
rf_out_err = NaN(ntrials,length(dims));
f1_out_err = NaN(ntrials,length(dims));
f2_out_err = NaN(ntrials,length(dims));

for trial = 1:ntrials
    
    fprintf('trial %d\n',trial)
    
    for i = 1:length(dims)
    
        d = dims(i);
        fprintf('d = %d\n',d)
        df = 10*d;   %degrees of freedom for inverse wishart
        nvartosample = ceil(d^(2/3));
        
        Mu = zeros(J,d);
        Sigma = zeros(d,d,J);
        
        for j = 1:J
        Mu(j,:) = mvnrnd(zeros(1,d),eye(d));
        Sigma(:,:,j) = iwishrnd(eye(d),df)*(df-d-1);
        end

        obj = gmdistribution(Mu,Sigma,p);
        [X,idx] = random(obj,n);
        Y = cellstr(num2str(Class(idx)));
        
        R = random_rotation(d);
        T = random_translation(n,d,-1,1);
        S = random_scaling(n,d,0,2);
        Sigma_outlier = 4*Sigma;
        X_rot = X*R;
        X_trans = X + T;
        X_scale = X.*S;
        outlier_model = gmdistribution(Mu,Sigma_outlier,p);
        [X_out,idx_out] = random(obj,0.2*n);
        X_out = cat(1,X,X_out);
        Y_out = cellstr(num2str(Class(idx_out)));
        Y_out = cat(1,Y,Y_out);

        rf = rpclassificationforest2(ntrees,X,Y,'nvartosample',nvartosample,'RandomForest',true);
        rf_err(trial,i) = oobpredict(rf,X,Y,'last');
        clear rf

        f1 = rpclassificationforest2(ntrees,X,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        f1_err(trial,i) = oobpredict(f1,X,Y,'last');
        clear f1
        
        f2 = rpclassificationforest2(ntrees,X,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new','Robust',true);
        f2_err(trial,i) = oobpredict(f2,X,Y,'last');
        clear f2

        rf_rot = rpclassificationforest2(ntrees,X_rot,Y,'nvartosample',nvartosample,'RandomForest',true);
        rf_rot_err(trial,i) = oobpredict(rf_rot,X_rot,Y,'last');
        clear rf_rot

        f1_rot = rpclassificationforest2(ntrees,X_rot,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        f1_rot_err(trial,i) = oobpredict(f1_rot,X_rot,Y,'last');
        clear f1_rot
        
        f2_rot = rpclassificationforest2(ntrees,X_rot,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new','Robust',true);
        f2_rot_err(trial,i) = oobpredict(f2_rot,X_rot,Y,'last');
        clear f2_rot

        rf_trans = rpclassificationforest2(ntrees,X_trans,Y,'nvartosample',nvartosample,'RandomForest',true);
        rf_trans_err(trial,i) = oobpredict(rf_trans,X_trans,Y,'last');
        clear rf_trans

        f1_trans = rpclassificationforest2(ntrees,X_trans,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        f1_trans_err(trial,i) = oobpredict(f1_trans,X_trans,Y,'last');
        clear f1_trans
        
        f2_trans = rpclassificationforest2(ntrees,X_trans,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new','Robust',true);
        f2_trans_err(trial,i) = oobpredict(f2_trans,X_trans,Y,'last');
        clear f2_trans

        rf_scale = rpclassificationforest2(ntrees,X_scale,Y,'nvartosample',nvartosample,'RandomForest',true);
        rf_scale_err(trial,i) = oobpredict(rf_scale,X_scale,Y,'last');
        clear rf_scale

        f1_scale = rpclassificationforest2(ntrees,X_scale,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        f1_scale_err(trial,i) = oobpredict(f1_scale,X_scale,Y,'last');
        clear f1_scale
        
        f2_scale = rpclassificationforest2(ntrees,X_scale,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new','Robust',true);
        f2_scale_err(trial,i) = oobpredict(f2_scale,X_scale,Y,'last');
        clear f2_scale
        
        rf_out = rpclassificationforest2(ntrees,X_out,Y_out,'nvartosample',nvartosample,'RandomForest',true);
        rf_out_err(trial,i) = oobpredict(rf_out,X_out,Y_out,'last');
        clear rf_out

        f1_out = rpclassificationforest2(ntrees,X_out,Y_out,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        f1_out_err(trial,i) = oobpredict(f1_out,X_out,Y_out,'last');
        clear f1_out
        
        f2_out = rpclassificationforest2(ntrees,X_out,Y_out,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new','Robust',true);
        f2_out_err(trial,i) = oobpredict(f2_out,X_out,Y_out,'last');
        clear f2_out
    end
end

save('Invariance_Multimodal.mat','R','rf_err','f1_err','f2_err','rf_rot_err',...
    'f1_rot_err','f2_rot_err','rf_trans_err','f1_trans_err','f2_trans_err',...
    'rf_scale_err','f1_scale_err','f2_scale_err','rf_out_err','f1_out_err','f2_out_err')

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
sem_rf_out = std(rf_out_err)/sqrt(ntrials);
sem_f1_out = std(f1_out_err)/sqrt(ntrials);
sem_f2_out = std(f2_out_err)/sqrt(ntrials);

Ynames = {'mean_rf_err' 'mean_f1_err' 'mean_f2_err'};
Enames = {'sem_rf' 'sem_f1' 'sem_f2'};
lspec = {'-bo','-rx','-gd'};
facespec = {'b','r','g'};
subplot(2,3,1)
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
    hold on
end
set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel(sprintf('OOB Error for %d Trees',ntrees))
legend('Random Forest','Sparse Randomer Forest w/ Mean Diff','Robust Sparse Randomer Forest w/ Mean Diff')
title('Multimodal')

Ynames = {'mean_rf_rot_err' 'mean_f1_rot_err' 'mean_f2_rot_err'};
Enames = {'sem_rf_rot' 'sem_f1_rot' 'sem_f2_rot'};
subplot(2,3,2)
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
    hold on
end
set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel(sprintf('OOB Error for %d Trees',ntrees))
legend('Random Forest','Sparse Randomer Forest w/ Mean Diff','Robust Sparse Randomer Forest w/ Mean Diff')
title('Rotated')

Ynames = {'mean_rf_trans_err' 'mean_f1_trans_err' 'mean_f2_trans_err'};
Enames = {'sem_rf_trans' 'sem_f1_trans' 'sem_f2_trans'};
subplot(2,3,3)
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
    hold on
end
set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel(sprintf('OOB Error for %d Trees',ntrees))
legend('Random Forest','Sparse Randomer Forest w/ Mean Diff','Robust Sparse Randomer Forest w/ Mean Diff')
title('Translated')

Ynames = {'mean_rf_scale_err' 'mean_f1_scale_err' 'mean_f2_scale_err'};
Enames = {'sem_rf_scale' 'sem_f1_scale' 'sem_f2_scale'};
subplot(2,3,4)
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
    hold on
end
set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel(sprintf('OOB Error for %d Trees',ntrees))
legend('Random Forest','Sparse Randomer Forest w/ Mean Diff','Robust Sparse Randomer Forest w/ Mean Diff')
title('Scaled')

Ynames = {'mean_rf_out_err' 'mean_f1_out_err' 'mean_f2_out_err'};
Enames = {'sem_rf_out' 'sem_f1_out' 'sem_f2_out'};
subplot(2,3,5)
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
    hold on
end
set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel(sprintf('OOB Error for %d Trees',ntrees))
legend('Random Forest','Sparse Randomer Forest w/ Mean Diff','Robust Sparse Randomer Forest w/ Mean Diff')
title('Outliers')

filename = 'Invariance_Multimodal';
save_fig(gcf,filename)