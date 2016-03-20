close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

Color = [0 1 1;0 1 0;1 0 1;1 0 0;0 0 0];

load Sparse_parity_data
ndims = length(dims);
ntrees = 500;
nmixs = 2:6;
NWorkers = 2;

for i = 3:3

    d = dims(i);
    
    if d <= 5
        mtrys = 1:d;
    else
        mtrys = ceil(d.^[0 1 2]);
    end
    
    if d >= 6
        nmixs = 2:6;
    else
        nmixs = 2:d;
    end 
    
    for trial = 1:1

        fprintf('trial %d\n',trial)

        for j = 1:length(mtrys)
            
            mtry = mtrys(j);

            fprintf('mtry = %d\n',mtry)
            
            poolobj = gcp('nocreate');
            if isempty(poolobj)
                parpool('local',NWorkers);
            end
            
            tic;
            rerf = rpclassificationforest(ntrees,X{i}(:,:,trial),...
                Y{i}(:,trial),'sparsemethod','sparse',...
                'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            
            nVarsCombined = NaN(1,rerf.nTrees);

            for t = 1:rerf.nTrees
                InternalNodes = rerf.Tree{t}.var~=0;
                nInternalNodes = sum(InternalNodes);
                nVarsCombined(t) = sum(cellfun(@nnz,rerf.Tree{t}.rpm(InternalNodes)))/...
                    sum(InternalNodes);
            end

            [~,err] = oobpredict2(rerf,X{i}(:,:,trial),Y{i}(:,trial));

            plot(nVarsCombined,err,'.','Color',Color(j,:))
            hold on
        end
    end
end

xlabel('Average Projection Density')
ylabel('OOB Error')
legend(['mtry = ',num2str(mtrys(1))],['mtry = ',num2str(mtrys(2))],['mtry = ',num2str(mtrys(3))])
title('Sparse Parity (n = 1000, p = 10)')
save_fig(gcf,[rerfPath 'RandomerForest/Figures/Lhat_vs_avg_density'])