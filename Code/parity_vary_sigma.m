%Parity

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
d =100;
ntrials = 10;
ntrees = 1000;
Sigmas = logspace(-1,1,5);
cumberr = NaN(ntrials,length(Sigmas));
cumf1err = NaN(ntrials,length(Sigmas));
cumf2err = NaN(ntrials,length(Sigmas));
cumf3err = NaN(ntrials,length(Sigmas));
embeddims = ceil(d^(2/3));

for trial = 1:ntrials
    trial
    nones = randi(d,n,1);
    Y = mod(nones,2);
    Ystr = cellstr(num2str(Y));
    Mu = sparse(n,d);

    for j = 1:n
            onesidx = randsample(1:d,nones(j),false);
            Mu(j,onesidx) = 1;
    end


    for i = 1:length(Sigmas)
        X = sparse(n,d);
        Sigma = Sigmas(i)*ones(1,d);
        for j = 1:n
            X(j,:) = mvnrnd(Mu(j,:),Sigma);
        end
        b = TreeBagger(ntrees,full(X),Y,'oobpred','on','NVarToSample',embeddims);
        berr = oobError(b);
        cumberr(trial,i) = berr(end);
        f1 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',embeddims,'mdiff','off','sparsemethod','new');
        f1err = oobpredict(f1,X,Ystr,'every');
        cumf1err(trial,i) = f1err(end);
        f2 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',embeddims,'mdiff','off','sparsemethod','old');
        f2err = oobpredict(f2,X,Ystr,'every');
        cumf2err(trial,i) = f2err(end);
        f3 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',embeddims,'mdiff','on','sparsemethod','new');
        f3err = oobpredict(f3,X,Ystr,'every');
        cumf3err(trial,i) = f3err(end);        
    end
    %plot(Sigmas,cumberr(trial,:),'bo',Sigmas,cumf1err(trial,:),'rx',Sigmas,cumf2err(trial,:),'gd',Sigmas,cumf3err(trial,:),'ks')
end
%save('xor.mat','cumberr','cumf1err','cumf2err','cumf3err')
bsem = std(cumberr)/sqrt(ntrials);
f1sem = std(cumf1err)/sqrt(ntrials);
f2sem  = std(cumf2err)/sqrt(ntrials);
f3sem = std(cumf3err)/sqrt(ntrials);
cumberr = mean(cumberr);
cumf1err = mean(cumf1err);
cumf2err = mean(cumf2err);
cumf3err = mean(cumf3err);
Ynames = {'cumberr' 'cumf1err' 'cumf2err' 'cumf3err'};
Enames = {'bsem' 'f1sem' 'f2sem' 'f3sem'};
lspec = {'-bo','-rx','-gd','-ks'};
facespec = {'b','r','g','k'};
hold on
for i = 1:length(Ynames)
    errorbar(Sigmas,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
end
set(gca,'XScale','log')
xlabel('variance')
ylabel('oob error')
legend('RandomForest','TylerForest+','TylerForest','TylerForest+meandiff')
%save_fig(gcf,sprintf('parity_ooberror_vs_sigma_n%d_d%d_ntrees%d',n,d,ntrees))