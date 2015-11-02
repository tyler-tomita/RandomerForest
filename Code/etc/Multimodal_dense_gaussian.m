clear
close all
clc

n = 100;    %# samples
dims = round(logspace(log10(2),3,5));
ntrials = 10;
ntrees = 1000;
nclasses = 4;
Classes = 1:nclasses;
ncomponents = 2;    %# of mixture components per class
J = nclasses*ncomponents; %total number of mixture components
p = ones(1,J)/J;    %mixture probabilities

%randomly allocate components evenly across all classes
Class = repmat(transpose(1:nclasses),1,ncomponents);
Class = Class(:);
Class = Class(randperm(length(Class)));

f_err = NaN(ntrials,length(dims));

for trial = 1:ntrials
    fprintf('trial %d\n',trial)
    
    for i = 1:length(dims)    
        fprintf('d = %d\n',dims(i))
        
        d = dims(i);
        nvartosample = ceil(d^(2/3));
        df = 10*d;   %degrees of freedom for inverse wishart
        
        Mu = zeros(J,d);
        Sigma = zeros(d,d,J);
        
        for j = 1:J
        Mu(j,:) = mvnrnd(zeros(1,d),eye(d));
        Sigma(:,:,j) = iwishrnd(eye(d),df)*(df-d-1);
        end

        obj = gmdistribution(Mu,Sigma,p);
        [X,idx] = random(obj,n);
        Y = cellstr(num2str(Class(idx)));
        
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
fname = 'Multimodal_dense_gaussian';
save_fig(gcf,fname)