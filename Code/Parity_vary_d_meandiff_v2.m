%Parity with varying # ambient dimensions

close all
clear
clc

n = 100;
dims = round(logspace(log10(2),3,5));
ntrials = 10;
ntrees = 1000;
ferr = NaN(ntrials,length(dims));

for trial = 1:ntrials
    fprintf('trial %0.0f\n',trial)

    for i = 1:length(dims)
        d = dims(i);
        fprintf('d = %0.0f\n',d)
        nvartosample = ceil(d^(2/3));
        X = sparse(n,d);
        Sigma = ones(1,d);
        nones = randi(d,n,1);
        Y = mod(nones,2);
        Ystr = cellstr(num2str(Y));
        Mu = sparse(n,d);
        for j = 1:n
            onesidx = randsample(1:d,nones(j),false);
            Mu(j,onesidx) = 1;
        end
        for j = 1:n
            X(j,:) = mvnrnd(Mu(j,:),Sigma);
        end

        f = rpclassificationforest2(ntrees,X,Ystr,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        ferr(trial,i) = oobpredict(f,X,Ystr,'last');
    end
end

save('Parity_vary_d_meandiff_v2.mat','ferr')
fsem = std(ferr)/sqrt(ntrials);
mean_ferr = mean(ferr);
errorbar(dims,mean_ferr,fsem);
set(gca,'XScale','log')
xlabel('# ambient dimensions')
ylabel('OOB Error')
xlabel('Number of Ambient Dimensions')
title('Parity')
legend('TylerForest+meandiff')
fname = sprintf('Parity/Parity_ooberror_vs_d_meandiff_v2_n%d_var%d_ntrees%d_ntrials%d',n,Sigma(1),ntrees,ntrials);
save_fig(gcf,fname)