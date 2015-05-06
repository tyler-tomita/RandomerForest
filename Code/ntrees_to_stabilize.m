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

n = 100;
ntrials = 10;
ntrees = 2000;
nclasses = 4;
Classes = 1:nclasses;
ncomponents = 2;    %# of mixture components per class
J = nclasses*ncomponents; %total number of mixture components
p = ones(1,J)/J;    %mixture probabilities

%randomly allocate components evenly across all classes
Class_mm = repmat(transpose(1:nclasses),1,ncomponents);
Class_mm = Class_mm(:);
Class_mm = Class_mm(randperm(length(Class_mm)));

cumrferr = NaN(ntrials,ntrees,3);
cumf1err = NaN(ntrials,ntrees,3);
cumf2err = NaN(ntrials,ntrees,3);
cumf3err = NaN(ntrials,ntrees,3);
Class_trunk = [0;1];
Title = {'Trunk' 'Parity' 'Multimodal'};

poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local',16,'IdleTimeout',120);
end

for trial = 1:ntrials
    fprintf('trial %d\n',trial)
    
    %Trunk
    d = 1000;
    d_idx = 1:d;
    nvartosample = ceil(d^(2/3));
    mu1 = 1./sqrt(d_idx);
    mu0 = -1*mu1;
    Mu = cat(1,mu0,mu1);
    Sigma = ones(1,d);
    obj = gmdistribution(Mu,Sigma);
    [X,idx] = random(obj,n);
    Y = cellstr(num2str(Class_trunk(idx)));
    
    rf = rpclassificationforest2(ntrees,X,Y,'nvartosample',nvartosample,'RandomForest',true);
    cumrferr(trial,:,1) = transpose(oobpredict(rf,X,Y,'every'));
    clear rf

    fprintf('Random Forest complete\n')
    
    f1 = rpclassificationforest2(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','old');
    cumf1err(trial,:,1) = transpose(oobpredict(f1,X,Y,'every'));
    clear f1
    
    fprintf('TylerForest complete\n')

    f2 = rpclassificationforest2(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','new');
    cumf2err(trial,:,1) = transpose(oobpredict(f2,X,Y,'every'));
    clear f2
    
    fprintf('TylerForest+ complete\n')
    
    f3 = rpclassificationforest2(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
    cumf3err(trial,:,1) = transpose(oobpredict(f3,X,Y,'every'));
    clear f3
    
    fprintf('TylerForest+meandiff complete\n')
    
    fprintf('Trunk complete\n\n')
    
    %Parity
    d = 10;
    X = sparse(n,d);
    Sigma = 1/32*speye(d);
    Mu = sparse(n,d);

    for j = 1:n
        Mu(j,:) = binornd(1,0.5,1,d);
        X(j,:) = mvnrnd(Mu(j,:),Sigma);
    end
    
    nones = sum(Mu,2);
    Y = cellstr(num2str(mod(nones,2)));
    
    rf = rpclassificationforest2(ntrees,X,Y,'nvartosample',nvartosample,'RandomForest',true);
    cumrferr(trial,:,2) = transpose(oobpredict(rf,X,Y,'every'));
    clear rf

    fprintf('Random Forest complete\n')
    
    f1 = rpclassificationforest2(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','old');
    cumf1err(trial,:,2) = transpose(oobpredict(f1,X,Y,'every'));
    clear f1

    fprintf('TylerForest complete\n')

    f2 = rpclassificationforest2(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','new');
    cumf2err(trial,:,2) = transpose(oobpredict(f2,X,Y,'every'));
    clear f2

    fprintf('TylerForest+ complete\n')
    
    f3 = rpclassificationforest2(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
    cumf3err(trial,:,2) = transpose(oobpredict(f3,X,Y,'every'));
    clear f3

    fprintf('TylerForest+meandiff complete\n')
    
    fprintf('Parity complete\n\n')
    
    %Multimodal
    d = 1000;
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
    Y = cellstr(num2str(Class_mm(idx)));
    
    rf = rpclassificationforest2(ntrees,X,Y,'nvartosample',nvartosample,'RandomForest',true);
    cumrferr(trial,:,3) = transpose(oobpredict(rf,X,Y,'every'));
    clear rf

    fprintf('Random Forest complete\n')
    
    f1 = rpclassificationforest2(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','old');
    cumf1err(trial,:,3) = transpose(oobpredict(f1,X,Y,'every'));
    clear f1

    fprintf('TylerForest complete\n')

    f2 = rpclassificationforest2(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','new');
    cumf2err(trial,:,3) = transpose(oobpredict(f2,X,Y,'every'));
    clear f2

    fprintf('TylerForest+ complete\n')
    
    f3 = rpclassificationforest2(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
    cumf3err(trial,:,3) = transpose(oobpredict(f3,X,Y,'every'));
    clear f3
    
    fprintf('TylerForest+meandiff complete\n')
    
    fprintf('Multimodal complete\n\n')
end

rfsem = std(cumrferr)/sqrt(ntrials);
f1sem = std(cumf1err)/sqrt(ntrials);
f2sem  = std(cumf2err)/sqrt(ntrials);
f3sem = std(cumf3err)/sqrt(ntrials);
cumrferr = mean(cumrferr);
cumf1err = mean(cumf1err);
cumf2err = mean(cumf2err);
cumf3err = mean(cumf3err);
Ynames = {'cumrferr(:,:,k)' 'cumf1err(:,:,k)' 'cumf2err(:,:,k)' 'cumf3err(:,:,k)'};
%Enames = {'rfsem(:,:,k)' 'f1sem(:,:,k)' 'f2sem(:,:,k)' 'f3sem(:,:,k)'};
lspec = {'-b','-r','-g','-k'};
for k = 1:length(Title)
    subplot(1,3,k)
    for i = 1:length(Ynames)
        plot(1:ntrees,eval(Ynames{i}),lspec{i});
        hold on
    end
    xlabel('# of Trees')
    ylabel('oob error')
    title(Title{k})
    legend('RandomForest','TylerForest','TylerForest+','TylerForest+meandiff')
end
fname = sprintf('ntrees_to_stabilize');
save_fig(gcf,fname)