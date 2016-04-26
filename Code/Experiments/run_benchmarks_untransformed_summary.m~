%% Compute oob error, ave tree strength, tree bias, and tree variance

close all
clear
clc

InPath = '../Results/Raw/Untransformed/';
contents = dir([InPath,'*.mat']);
OutPath = '../Results/Summary/Untransformed/';

load('../Data/Benchmark_data.mat','Y','n','d')

min_d = 1;
max_d = 100;
min_n = 1;
max_n = 50000;

rmidx = [];
for i = 1:length(n)
    if ~(n(i) >= min_n && n(i) < max_n && d(i) >= min_d && d(i) < max_d)
        rmidx = [rmidx i];
    end
end
n(rmidx) = [];
d(rmidx) = [];
Y(rmidx) = [];

for i = 1:length(contents)
    InFile = [InPath,contents(i).name];
    load(InFile,'Yhats','Time','ntrees','mtrys','nmixs')
    Classifiers = fieldnames(Yhats);
    mtrys_rf = mtrys(mtrys<=d(i));
    
    for c = 1:length(Classifiers)
        cl = Classifiers{c};
        if strcmp(cl,'rf') || strcmp(cl,'rf_rot')
            m = mtrys_rf;
        else
            m = mtrys;
        end
        
        MR.(cl) = [];
        S.(cl) = [];
        V.(cl) = [];
        
        if ndims(Yhats.(cl)) == 3
            for j = 1:length(m)
                MR.(cl)(j,1) = oob_error(Yhats.(cl)(:,:,j),Y{i},'last');
                S.(cl)(j,1) = misclassification_rate(Yhats.(cl)(:,:,j),...
                    Y{i},false);
                V.(cl)(j,1) = classifier_variance(Yhats.(cl)(:,:,j));
            end
            NMIX.(cl) = [];
        else
            for j = 1:length(m)
                for k = 1:length(nmixs)
                    MR.(cl)(j,k) = oob_error(Yhats.(cl)(:,:,j,k),Y{i},'last');
                    S.(cl)(j,k) = misclassification_rate(Yhats.(cl)(:,:,j,k),...
                        Y{i},false);
                    V.(cl)(j,k) = classifier_variance(Yhats.(cl)(:,:,j,k));
                end
            end
            NMIX.(cl) = nmixs;
        end
        B.(cl) = S.(cl) - V.(cl);
        MTRY.(cl) = m;
    end
    Summary = struct('MR',MR,'S',S,'V',V,'B',B,'MTRY',MTRY,'NMIX',NMIX);
    save([OutPath,contents(i).name(1:end-4),'_untransformed_summary.mat'],...
        'Summary')
end