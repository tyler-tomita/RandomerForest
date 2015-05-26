%Parity with varying # ambient dimensions

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
dims = [2 3 4 5 10];
ntrials = 10;
ntrees = 1000;
NWorkers = 16;
cumrferr = NaN(ntrials,length(dims));
cumf1err = NaN(ntrials,length(dims));
cumf2err = NaN(ntrials,length(dims));
cumf3err = NaN(ntrials,length(dims));
cumf4err = NaN(ntrials,length(dims));
trf = NaN(ntrials,length(dims));
tf1 = NaN(ntrials,length(dims));
tf2 = NaN(ntrials,length(dims));
tf3 = NaN(ntrials,length(dims));
tf4 = NaN(ntrials,length(dims));
for trial = 1:ntrials
    fprintf('trial %d\n',trial)

    for i = 1:length(dims)
        d = dims(i);
        fprintf('d = %d\n',d)
        nvartosample = ceil(d^(2/3));
        X = zeros(n,d);
        %Sigma = 1/8*ones(1,d);
        Sigma = 1/32*ones(1,d);
        %nones = randi(d+1,n,1)-1;
        %Y = mod(nones,2);
        %Ystr = cellstr(num2str(Y));
        Mu = sparse(n,d);
        for j = 1:n
            %onesidx = randsample(1:d,nones(j),false);
            %Mu(j,onesidx) = 1;
            Mu(j,:) = binornd(1,0.5,1,d);
            X(j,:) = mvnrnd(Mu(j,:),Sigma);
        end
        
        nones = sum(Mu,2);
        Y = cellstr(num2str(mod(nones,2)));

        tic
        rf = rpclassificationforest(ntrees,X,Y,'nvartosample',nvartosample,'RandomForest',true,'NWorkers',NWorkers);
        trf(trial,i) = toc;
        cumrferr(trial,i) = oobpredict(rf,X,Y,'last');
        
        fprintf('Random Forest complete\n')
        
        tic
        f1 = rpclassificationforest(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','dense','NWorkers',NWorkers);
        tf1(trial,i) = toc;
        cumf1err(trial,i) = oobpredict(f1,X,Y,'last');
        
        fprintf('TylerForest complete\n')
        
        tic
        f2 = rpclassificationforest(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','sparse','NWorkers',NWorkers);
        tf2(trial,i) = toc;
        cumf2err(trial,i) = oobpredict(f2,X,Y,'last');
        
        fprintf('TylerForest+ complete\n')
        
        tic
        f3 = rpclassificationforest(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','on','sparsemethod','sparse','NWorkers',NWorkers);
        tf3(trial,i) = toc;
        cumf3err(trial,i) = oobpredict(f3,X,Y,'last');
        
        fprintf('TylerForest+meandiff complete\n')
        
        tic
        f4 = rpclassificationforest(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','on','sparsemethod','sparse','Robust',true,'NWorkers',NWorkers);
        tf4(trial,i) = toc;
        cumf4err(trial,i) = oobpredict(f4,X,Y,'last');
        
        fprintf('Robust complete\n')
    end
end

save('Parity_vary_d.mat','cumrferr','cumf1err','cumf2err','cumf3err','cumf4err','trf','tf1','tf2','tf3','tf4')

rfsem = std(cumrferr)/sqrt(ntrials);
f1sem = std(cumf1err)/sqrt(ntrials);
f2sem  = std(cumf2err)/sqrt(ntrials);
f3sem  = std(cumf3err)/sqrt(ntrials);
f4sem  = std(cumf4err)/sqrt(ntrials);
cumrferr = mean(cumrferr);
cumf1err = mean(cumf1err);
cumf2err = mean(cumf2err);
cumf3err = mean(cumf3err);
cumf4err = mean(cumf4err);
Ynames = {'cumrferr' 'cumf1err' 'cumf2err' 'cumf3err' 'cumf4err'};
Enames = {'rfsem' 'f1sem' 'f2sem' 'f3sem' 'f4sem'};
lspec = {'-bo','-rx','-gd','-ks' '-m.'};
facespec = {'b','r','g','k','m'};
hold on
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
end
%set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel(sprintf('OOB Error for %d Trees',ntrees))
legend('RandomForest','TylerForest','TylerForest+','TylerForest+meandiff','Robust TylerForest+meandiff')
fname = sprintf('Parity_ooberror_vs_d_n%d_ntrees%d_ntrials%d_v6',n,ntrees,ntrials);
save_fig(gcf,fname)

rfsem = std(trf)/sqrt(ntrials);
f1sem = std(tf1)/sqrt(ntrials);
f2sem = std(tf2)/sqrt(ntrials);
f3sem = std(tf3)/sqrt(ntrials);
f4sem = std(tf4)/sqrt(ntrials);
trf = mean(trf);
tf1 = mean(tf1);
tf2 = mean(tf2);
tf3 = mean(tf3);
tf4 = mean(tf4);
Ynames = {'trf' 'tf1' 'tf2' 'tf3' 'tf4'};

figure(2)
hold on
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
end
%set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel('Wall Time (sec)')
legend('RandomForest','TylerForest','TylerForest+','TylerForest+meandiff','Robust TylerForest+meandiff')
fname = sprintf('Parity_time_vs_d_n%d_ntrees%d_ntrials%d_v6',n,ntrees,ntrials);
save_fig(gcf,fname)
