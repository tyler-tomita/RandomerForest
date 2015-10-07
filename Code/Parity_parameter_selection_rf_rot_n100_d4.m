close all
clear
clc

n = 100;
d =4;
ntrees = 1000;
ntrials = 20;
NWorkers = 2;
Colors = linspecer(4,'sequential');
mtrys = [1 2 3 4];
err_rf_rot = zeros(ntrees,length(mtrys),ntrials);

for trial = 1:ntrials
    
    fprintf('trial %d\n',trial)
    
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
    Ystr = cellstr(num2str(mod(nones,2)));    
    
    i = 1;
    
    for mtry = mtrys
        
        fprintf('mtry = %d\n',mtry)
        
        rf_rot = rpclassificationforest(ntrees,X,Ystr,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'rotate',true);
        err_rf_rot(:,i,trial) = oobpredict(rf_rot,X,Ystr,'every');
        
        i = i + 1;
    end
end

save(sprintf('Parity_parameter_selection_rf_rot_n%d_d%d.mat',n,d),'err_rf_rot')

mean_err_rf_rot = mean(err_rf_rot,3);

i = 1;

for mtry = mtrys
    plot(mean_err_rf_rot(:,i),'Color',Colors(i,:))
    hold on
    i = i + 1;
end

legend('rf rotate (mtry = 1)','rf rotate (mtry = 2)','rf rotate (mtry = 3)','rf rotate (mtry = 4)')

xlabel('# trees')
ylabel('oob error')
title(sprintf('Parity (n=%d, d=%d, ntrials=%d)',n,d,ntrials))

save_fig(gcf,sprintf('Parity_parameter_selection_rf_rot_n%d_d%d',n,d))
