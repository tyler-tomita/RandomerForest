close all
clear
clc

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

ferr = NaN(ntrials,length(dims));

for trial = 1:ntrials
    fprintf('trial %d\n',trial)
    
    for i = 1:length(dims)    
        fprintf('d = %d\n',dims(i))
        
        d = dims(i);
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
        Y = Class(idx);
        Ystr = cellstr(num2str(Y));

        f = rpclassificationforest2(ntrees,X,Ystr,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        ferr(trial,i) = oobpredict(f,X,Ystr,'last');
    end
end

save('Correlated_gaussians/Multimodal_vary_d_meandiff_v2.mat','ferr')
fsem = std(ferr)/sqrt(ntrials);
mean_ferr = mean(ferr);
errorbar(dims,mean_ferr,fsem);
set(gca,'XScale','log')
xlabel('# Ambient Dimensions')
ylabel(sprintf('OOB Error for %d Trees',ntrees))
legend('TylerForest+meandiff')
fname = sprintf('Correlated_gaussians/Multimodal_ooberror_vs_d_meandiff_v2_n%d_var%d_ntrees%d_ntrials%d',n,1,ntrees,ntrials);
save_fig(gcf,fname)
