function [Lhat,Time,n,d,sp,tr] = Benchmark(FileName,ntrees,NWorkers,ntrials)

    X = dlmread(FileName,'\t',1,1);
    Y = cellstr(num2str(X(:,end)));
    X(:,end) = [];
    [n,d] = size(X);
    tr = trace(cov(X));
    if d <= 5
        mtrys = 1:d;
    else
        mtrys = ceil(d.^[0 1/4 1/2 3/4 1]);
    end
    
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool('local',NWorkers,'IdleTimeout',360);
    end
    
    for trial = 1:ntrials
        
        i = 1;

        for mtry = mtrys

            tic
            rf = rpclassificationforest(ntrees,X,Y,'nvartosample',mtry,'RandomForest',true,'NWorkers',NWorkers);
            Time.rf(trial,i) = toc;
            Lhat.rf(:,i,trial) = oobpredict(rf,X,Y,'every');
            sp.rf(trial,i) = sum(rf.NumVars);
            clear rf

            tic
            rerf = rpclassificationforest(ntrees,X,Y,'nvartosample',mtry,'sparsemethod','sparse','NWorkers',NWorkers);
            Time.rerf(trial,i) = toc;
            Lhat.rerf(:,i,trial) = oobpredict(rerf,X,Y,'every');
            sp.rerf(trial,i) = sum(rerf.NumVars);
            clear rerf

            tic
            rf_rot = rpclassificationforest(ntrees,X,Y,'nvartosample',mtry,'RandomForest',true,'NWorkers',NWorkers,'rotate',true);
            Time.rf_rot(trial,i) = toc;
            Lhat.rf_rot(:,i,trial) = oobpredict(rf_rot,X,Y,'every');
            sp.rf_rot(trial,i) = sum(rf_rot.NumVars);
            clear rf_rot

            tic
            rerf_rot = rpclassificationforest(ntrees,X,Y,'nvartosample',mtry,'sparsemethod','sparse','NWorkers',NWorkers,'rotate',true);
            Time.rerf_rot(trial,i) = toc;
            Lhat.rerf_rot(:,i,trial) = oobpredict(rerf_rot,X,Y,'every');
            sp.rerf_rot(trial,i) = sum(rerf_rot.NumVars);
            clear rerf_rot

            i = i + 1;
        end
    end
end
