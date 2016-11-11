%% Plot Sparse Parity Posterior Heat Maps 
% Trains RF, F-RC, and Rotation RF on sparse parity and plots
% posterior heat maps for each classifier

%% Initialize parameters
close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Sparse_parity_data
load Random_matrix_adjustment_factor

Classifiers = {'rerf' 'rerfr' 'rerf2' 'rerf2r'};

Transformations = {'Untransformed' 'Affine' 'Outlier'};

ntrials = 10;

xmin = -1;
xmax = 1;
ymin = xmin;
ymax = xmax;
npoints = 50;
[xgv,ygv] = meshgrid(linspace(xmin,xmax,npoints),linspace(ymin,ymax,npoints));
Xpost = xgv(:);
Ypost = ygv(:);


for i = 3:3
    p = dims(i);
    fprintf('p = %d\n',p)
    
    
    Zpost = -0.5*ones(npoints^2,p-2);
    
    if p <= 10
        dx = round(p^3.5);
    elseif p > 10 && p <= 100
        dx = p^2;
    elseif p > 100 && p <= 1000
        dx = round(p^1.5);
    else
        dx = 10*p;
    end
    
    if dx <= 5
        mtrys_rerf2 = 1:dx;
    else
        mtrys_rerf2 = ceil(dx.^[1/4 1/2 3/4 1]);
    end
      
    if p <= 5
        mtrys = [1:p ceil(p.^[1.5 2])];
    elseif p > 5 && p <= 20
        mtrys = ceil(p.^[1/4 1/2 3/4 1 1.5 2]);
    else
        mtrys = [ceil(p.^[1/4 1/2 3/4 1]) 5*p 10*p];
    end
    mtrys_rf = mtrys(mtrys<=p);

    for c = 1:length(Classifiers)
        fprintf('%s start\n',Classifiers{c})
        Params{i}.(Classifiers{c}).nTrees = 500;
        Params{i}.(Classifiers{c}).Stratified = true;
        Params{i}.(Classifiers{c}).NWorkers = 2;
        if strcmp(Classifiers{c},'rfr') || strcmp(Classifiers{c},...
                'rerfr') || strcmp(Classifiers{c},'frcr') || ...
                strcmp(Classifiers{c},'rr_rfr') || strcmp(Classifiers{c},'rerf2r')
            Params{i}.(Classifiers{c}).Rescale = 'rank';
%         elseif strcmp(Classifiers{c},'rfn') || strcmp(Classifiers{c},...
%                 'rerfn') || strcmp(Classifiers{c},'frcn') || ...
%                 strcmp(Classifiers{c},'rr_rfn')
%             Params{i}.(Classifiers{c}).Rescale = 'normalize';
%         elseif strcmp(Classifiers{c},'rfz') || strcmp(Classifiers{c},...
%                 'rerfz') || strcmp(Classifiers{c},'frcz') || ...
%                 strcmp(Classifiers{c},'rr_rfz')
%             Params{i}.(Classifiers{c}).Rescale = 'zscore';
        else
            Params{i}.(Classifiers{c}).Rescale = 'off';
        end
        if strcmp(Classifiers{c},'rerfd')
            Params{i}.(Classifiers{c}).mdiff = 'node';
        else
            Params{i}.(Classifiers{c}).mdiff = 'off';
        end
        if strcmp(Classifiers{c},'rf') || strcmp(Classifiers{c},'rfr')...
                || strcmp(Classifiers{c},'rfn') || strcmp(Classifiers{c},'rfz') || ...
                strcmp(Classifiers{c},'rr_rf') || strcmp(Classifiers{c},'rr_rfr') || ...
                strcmp(Classifiers{c},'rr_rfn') || strcmp(Classifiers{c},'rr_rfz')
            Params{i}.(Classifiers{c}).ForestMethod = 'rf';
            Params{i}.(Classifiers{c}).d = mtrys_rf;
        elseif strcmp(Classifiers{c},'rerf') || strcmp(Classifiers{c},'rerfr')...
                || strcmp(Classifiers{c},'rerfn') || strcmp(Classifiers{c},'rerfz') || ...
                strcmp(Classifiers{c},'rerfd')
            Params{i}.(Classifiers{c}).ForestMethod = 'sparse-binary';
            Params{i}.(Classifiers{c}).d = mtrys;
            for j = 1:length(Params{i}.(Classifiers{c}).d)
                Params{i}.(Classifiers{c}).dprime(j) = ...
                    ceil(Params{i}.(Classifiers{c}).d(j)^(1/interp1(dims,...
                    slope,p)));
            end
        elseif strcmp(Classifiers{c},'frc') || strcmp(Classifiers{c},'frcr') || ...
                strcmp(Classifiers{c},'frcn') || strcmp(Classifiers{c},'frcz')
            Params{i}.(Classifiers{c}).ForestMethod = 'frc';
            Params{i}.(Classifiers{c}).d = mtrys;
            Params{i}.(Classifiers{c}).nmix = 2;
        elseif strcmp(Classifiers{c},'rerf2') || strcmp(Classifiers{c},'rerf2r')
            Params{i}.(Classifiers{c}).ForestMethod = 'rerf2';
            Params{i}.(Classifiers{c}).d = mtrys_rerf2;
            Params{i}.(Classifiers{c}).dx = dx;
        end
        if strcmp(Classifiers{c},'rr_rf') || strcmp(Classifiers{c},'rr_rfr') || ...
                strcmp(Classifiers{c},'rr_rfn') || strcmp(Classifiers{c},'rr_rfz')
            Params{i}.(Classifiers{c}).Rotate = true;
        end
        
        for t = 1:length(Transformations)
            fprintf('%s\n',Transformations{t})

            OOBError{i}.(Classifiers{c}).(Transformations{t}) = NaN(ntrials,length(Params{i}.(Classifiers{c}).d));
            OOBAUC{i}.(Classifiers{c}).(Transformations{t}) = NaN(ntrials,length(Params{i}.(Classifiers{c}).d));
            TrainTime{i}.(Classifiers{c}).(Transformations{t}) = NaN(ntrials,length(Params{i}.(Classifiers{c}).d));

            for trial = 1:ntrials
                fprintf('Trial %d\n',trial)

                % train classifier
                poolobj = gcp('nocreate');
                if isempty(poolobj)
                    parpool('local',Params{i}.(Classifiers{c}).NWorkers,...
                        'IdleTimeout',360);
                end

                tic;
                [Forest,~,TrainTime{i}.(Classifiers{c}).(Transformations{t})(trial,:)] = ...
                    RerF_train(Xtrain(i).(Transformations{t})(:,:,trial),...
                    Ytrain(i).(Transformations{t})(:,trial),Params{i}.(Classifiers{c}));

                % select best hyperparameter

                for j = 1:length(Params{i}.(Classifiers{c}).d)
                    Scores = rerf_oob_classprob(Forest{j},...
                        Xtrain(i).(Transformations{t})(:,:,trial),'last');
                    Predictions = predict_class(Scores,Forest{j}.classname);
                    OOBError{i}.(Classifiers{c}).(Transformations{t})(trial,j) = ...
                        misclassification_rate(Predictions,Ytrain(i).(Transformations{t})(:,trial),...
                        false);
                    if size(Scores,2) > 2
                        Yb = binarize_labels(Ytrain(i).(Transformations{t})(:,trial),Forest{j}.classname);
                        [~,~,~,OOBAUC{i}.(Classifiers{c}).(Transformations{t})(trial,j)] = ...
                            perfcurve(Yb(:),Scores(:),'1');
                    else
                        [~,~,~,OOBAUC{i}.(Classifiers{c}).(Transformations{t})(trial,j)] = ...
                            perfcurve(Ytrain(i).(Transformations{t})(:,trial),Scores(:,2),'1');
                    end
                end
                BestIdx = hp_optimize(OOBError{i}.(Classifiers{c}).(Transformations{t})(trial,:),...
                    OOBAUC{i}.(Classifiers{c}).(Transformations{t})(trial,:));
                if length(BestIdx)>1
                    BestIdx = BestIdx(end);
                end

                x = [Xpost,Ypost,Zpost];
                if strcmp(Transformations{t},'Affine')
                    x = (x*R{i}(:,:,trial)).*repmat(S2{i}(1,:,trial),length(Xpost),1);
                end
                if strcmp(Forest{BestIdx}.Rescale,'off')
                    Phats{i}.(Classifiers{c}).(Transformations{t})(:,:,trial) ...
                        = rerf_classprob(Forest{BestIdx},...
                        x,'last');
                else
                    Phats{i}.(Classifiers{c}).(Transformations{t})(:,:,trial) ...
                        = rerf_classprob(Forest{BestIdx},...
                        x,'last',...
                        Xtrain(i).(Transformations{t})(:,:,trial));
                end
                
                save([rerfPath 'RandomerForest/Results/Sparse_parity_uniform_transformations_posteriors_RerF.mat'],'dims',...
                    'Phats','Xpost','Ypost','Zpost')
            end
        end
        fprintf('%s complete\n',Classifiers{c})
    end   
end