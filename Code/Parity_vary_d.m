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
cumrferr = NaN(ntrials,length(dims));
cumf1err = NaN(ntrials,length(dims));
cumf2err = NaN(ntrials,length(dims));
cumf3err = NaN(ntrials,length(dims));
trf = NaN(ntrials,length(dims));
tf1 = NaN(ntrials,length(dims));
tf2 = NaN(ntrials,length(dims));
tf3 = NaN(ntrials,length(dims));

parpool;
for trial = 1:ntrials
    fprintf('trial %d\n',trial)

    for i = 1:length(dims)
        d = dims(i);
        fprintf('d = %d\n',d)
        nvartosample = ceil(d^(2/3));
        X = sparse(n,d);
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
        Y = mod(nones,2);
        Ystr = cellstr(num2str(Y));

        tic
        rf = rpclassificationforest(ntrees,X,Ystr,'nvartosample',nvartosample,'RandomForest',true);
        trf(trial,i) = toc;
        cumrferr(trial,i) = oobpredict(rf,X,Ystr,'last');
        clear rf


        fprintf('Random Forest complete\n')

        tic
        f1 = rpclassificationforest(ntrees,X,Ystr,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','old');
        tf1(trial,i) = toc;
        cumf1err(trial,i) = oobpredict(f1,X,Ystr,'last');
        clear f1

        fprintf('TylerForest complete\n')

        tic
        f2 = rpclassificationforest(ntrees,X,Ystr,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','new');
        tf2(trial,i) = toc;
        cumf2err(trial,i) = oobpredict(f2,X,Ystr,'last');
        clear f2

        fprintf('TylerForest+ complete\n')

        tic
        f3 = rpclassificationforest(ntrees,X,Ystr,'s',3,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        tf3(trial,i) = toc;
        cumf3err(trial,i) = oobpredict(f3,X,Ystr,'last');
        clear f3

        fprintf('TylerForest+meandiff complete\n')

    end
end

save('Parity_vary_d.mat','cumrferr','cumf1err','cumf2err','cumf3err','trf','tf1','tf2','tf3')

rfsem = std(cumrferr)/sqrt(ntrials);
f1sem = std(cumf1err)/sqrt(ntrials);
f2sem  = std(cumf2err)/sqrt(ntrials);
f3sem  = std(cumf3err)/sqrt(ntrials);
cumrferr = mean(cumrferr);
cumf1err = mean(cumf1err);
cumf2err = mean(cumf2err);
cumf3err = mean(cumf3err);
Ynames = {'cumrferr' 'cumf1err' 'cumf2err' 'cumf3err'};
Enames = {'rfsem' 'f1sem' 'f2sem' 'f3sem'};
lspec = {'-bo','-rx','-gd','-ks'};
facespec = {'b','r','g','k'};
hold on
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
end
%set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel(sprintf('OOB Error for %d Trees',ntrees))
legend('RandomForest','TylerForest','TylerForest+','TylerForest+meandiff')
fname = sprintf('Parity_ooberror_vs_d_n%d_ntrees%d_ntrials%d_v5',n,ntrees,ntrials);
save_fig(gcf,fname)

rfsem = std(trf)/sqrt(ntrials);
f1sem = std(tf1)/sqrt(ntrials);
f2sem = std(tf2)/sqrt(ntrials);
f3sem = std(tf3)/sqrt(ntrials);
trf = mean(trf);
tf1 = mean(tf1);
tf2 = mean(tf2);
tf3 = mean(tf3);
Ynames = {'trf' 'tf1' 'tf2' 'tf3'};

figure(2)
hold on
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
end
%set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel('Wall Time (sec)')
legend('RandomForest','TylerForest','TylerForest+','TylerForest+meandiff')
fname = sprintf('Parity_time_vs_d_n%d_ntrees%d_ntrials%d_v5',n,ntrees,ntrials);
save_fig(gcf,fname)
