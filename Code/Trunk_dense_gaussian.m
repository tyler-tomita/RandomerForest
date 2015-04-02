clear
close all
clc

n = 100;
dims = round(logspace(log10(2),3,5));
ntrees = 1000;
ntrials = 10;
f_err = NaN(ntrials,length(dims));
Class = [0;1];

for trial = 1:ntrials
    fprintf('trial %d\n',trial)

    for i = 1:length(dims)
        d = dims(i);
        fprintf('dimensions = %d\n',dims(i))
        d_idx = 1:d;
        mu1 = 1./sqrt(d_idx);
        mu0 = -1*mu1;
        Mu = cat(1,mu0,mu1);
        Sigma = 1*speye(d);
        obj = gmdistribution(Mu,Sigma);
        [X,idx] = random(obj,n);
        Y = cellstr(num2str(Class(idx)));
        f = rpclassificationforest2(ntrees,X,Y,'nvartosample',ceil(size(X,2)^(2/3)),'sparsemethod','dgaussian','s',3,'mdiff','on');
        f_err(trial,i) = oobpredict(f,X,Y,'last');
    end
end

mean_f_err = mean(f_err);
sem_f_err = std(f_err)/sqrt(ntrials);

errorbar(dims,mean_f_err,sem_f_err);
set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel('oob error')
%ylabel('Bayes error')
%save_fig(gcf,fname)