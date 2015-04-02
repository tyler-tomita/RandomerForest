clear
close all
clc

n = 100;
dims = round(logspace(log10(2),3,5));
ntrials = 10;
ntrees = 1000;
f_err = NaN(ntrials,length(dims));

for trial = 1:ntrials
    fprintf('trial %d\n',trial)

    for i = 1:length(dims)
        d = dims(i);
        fprintf('d = %d\n',d)
        nvartosample = ceil(d^(2/3));
        X = sparse(n,d);
        Sigma = 1/8*speye(d);
        %nones = randi(d,n,1);
        %Y = mod(nones,2);
        Mu = sparse(n,d);
        
        for j = 1:n
            %onesidx = randsample(1:d,nones(j),false);
            %Mu(j,onesidx) = 1;
            Mu(j,:) = binornd(1,0.5,1,d);
            X(j,:) = mvnrnd(Mu(j,:),Sigma);
        end
        nones = sum(Mu,2);
        Y = cellstr(num2str(mod(nones,2)));
        
        f = rpclassificationforest2(ntrees,X,Y,'nvartosample',nvartosample,'sparsemethod','dgaussian','s',3);
        f_err(trial,i) = oobpredict(f,X,Y,'last');

    end
end

mean_f_err = mean(f_err);
sem_f_err = std(f_err)/sqrt(ntrials);

errorbar(dims,mean_f_err,sem_f_err);
set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel('oob error')
fname = 'Parity_dense_gaussian';
save_fig(gcf,fname)