close all
clear
clc

n = 100;
d = 1000;
ntrees = 1000;
ntrials = 100;
embeddims = ceil([0.5 1 2 4 8]*sqrt(d));
cumberr = NaN(ntrials,length(embeddims));
cumf1err = NaN(ntrials,length(embeddims));
cumf2err = NaN(ntrials,length(embeddims));
cumf3err = NaN(ntrials,length(embeddims));
tb = NaN(ntrials,length(embeddims));
tf1 = NaN(ntrials,length(embeddims));
tf2 = NaN(ntrials,length(embeddims));
tf3 = NaN(ntrials,length(embeddims));

for trial = 1:ntrials
    fprintf('trial %d\n',trial)
    Mu0 = randn(1,d);
    L = randn(d);
    Sigma0 = L*L';
    X0 = mvnrnd(Mu0,Sigma0,round(n/2));
    Y0 = zeros(round(n/2),1);
    
    Mu1 = randn(1,d);
    L = randn(d);
    Sigma1 = L*L';
    X1 = mvnrnd(Mu1,Sigma1,round(n/2));
    Y1 = ones(round(n/2),1);
    
    X = cat(1,X0,X1);
    Y = cat(1,Y0,Y1);
    Ystr = cellstr(num2str(Y));
    
    for i = 1:length(embeddims)
        nvartosample = embeddims(i);
        fprintf('nvartosample = %d\n',nvartosample)
        poolobj = gcp('nocreate');
        if isempty(poolobj)
            parpool;
        end
        paroptions = statset('UseParallel',true);
        tic
        b = TreeBagger(ntrees,full(X),Y,'oobpred','on','NVarToSample',nvartosample,'Options',paroptions);
        tb(trial,i) = toc;
        berr = oobError(b);
        clear b
        cumberr(trial,i) = berr(end);
        

        fprintf('Random Forest complete\n')
        
        tic
        f1 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','old');
        tf1(trial,i) = toc;
        f1err = oobpredict(f1,X,Ystr,'every');
        clear f1
        cumf1err(trial,i) = f1err(end);
        
        fprintf('TylerForest complete\n')
        
        tic
        f2 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','new');
        tf2(trial,i) = toc;
        f2err = oobpredict(f2,X,Ystr,'every');
        clear f2
        cumf2err(trial,i) = f2err(end);
        
        fprintf('TylerForest+ complete\n')
        
        tic
        f3 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        tf3(trial,i) = toc;
        f3err = oobpredict(f3,X,Ystr,'every');
        clear f3
        cumf3err(trial,i) = f3err(end);
        
        fprintf('TylerForest+meandiff complete\n')
    end
end

save('Correlated_gaussian_2_class.mat','cumberr','cumf1err','cumf2err','cumf3err','tb','tf1','tf2','tf3')
bsem = std(cumberr)/sqrt(ntrials);
f1sem = std(cumf1err)/sqrt(ntrials);
f2sem  = std(cumf2err)/sqrt(ntrials);
f3sem  = std(cumf3err)/sqrt(ntrials);
cumberr = mean(cumberr);
cumf1err = mean(cumf1err);
cumf2err = mean(cumf2err);
cumf3err = mean(cumf3err);
Ynames = {'cumberr' 'cumf1err' 'cumf2err' 'cumf3err'};
Enames = {'bsem' 'f1sem' 'f2sem' 'f3sem'};
lspec = {'-bo','-rx','-gd','ks'};
facespec = {'b','r','g','k'};
hold on
for i = 1:length(Ynames)
    errorbar(embeddims,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
end
set(gca,'XScale','log')
xlabel('# of Embedded Dimensions')
ylabel(sprintf('OOB Error for %d Trees',ntrees))
legend('RandomForest','TylerForest','TylerForest+','TylerForest+meandiff')
fname = sprintf('correlated_gaussian_2_class_ooberror_vs_embeddims_n%d_d%d_ntrees%d_ntrials%d',size(X,1),d,ntrees,ntrials);
save_fig(gcf,fname)
