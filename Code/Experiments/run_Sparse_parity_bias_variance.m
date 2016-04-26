%% Computes estimates of bias, variance, and generalization error for different decision forest classifiers

%   Predictions.(classifier name) is an n.test-by-n.trainSet matrix of
%   predicted class labels for a classifier for each value of mtry
%
%   B.(classifier name) is a 1-by-ndims cell array of 1d arrays of bias
%   estimates for each value of mtry
%
%   V.(classifier name) is a 1-by-ndims cell array of 1d arrays of variance
%   estimates for each value of mtry
%
%   GE.(classifier name) is a 1-by-ndims cell array of 1d arrays of
%   estimates of generalization error for each value of mtry
%
%   Params.(classifier name) is a structure with fields RandomForest,
%   SparseMethod, mtry, nmix, ntrees, Stratified, NWorkers, Rotate

close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Sparse_parity_partitioned_data

ntrees = 500;
NWorkers = 16;
Stratified = true;
Classifiers = {'rf' 'rerf' 'rf_rot' 'frc2' 'frc3' 'frc4' 'frc5' 'frc6'};
ParameterNames = {'RandomForest' 'SparseMethod' 'mtry' 'nmix' 'ntrees' ...
    'Stratified' 'NWorkers' 'Rotate'};
RandomForestParams = [true false true false(1,5)];
SparseMethodParams = {'' 'sparse' '' 'frc' 'frc' 'frc' 'frc' 'frc'};
RotateParams = {false false true false(1,5)};

Params.rf.RandomForest = true;
Params.rf.SparseMethod = '';
Params.rf.Rotate = false;
Params.rf.nmix = [];

Params.rerf.RandomForest = false;
Params.rerf.SparseMethod = 'sparse';
Params.rerf.Rotate = false;
Params.rerf.nmix = [];

Params.rf_rot.RandomForest = true;
Params.rf_rot.SparseMethod = '';
Params.rf_rot.Rotate = true;
Params.rf_rot.nmix = [];

Params.frc2.RandomForest = false;
Params.frc2.SparseMethod = 'frc';
Params.frc2.Rotate = false;
Params.frc2.nmix = 2;

Params.frc3.RandomForest = false;
Params.frc3.SparseMethod = 'frc';
Params.frc3.Rotate = false;
Params.frc3.nmix = 3;

Params.frc4.RandomForest = false;
Params.frc4.SparseMethod = 'frc';
Params.frc4.Rotate = false;
Params.frc4.nmix = 4;

Params.frc5.RandomForest = false;
Params.frc5.SparseMethod = 'frc';
Params.frc5.Rotate = false;
Params.frc5.nmix = 5;

Params.frc6.RandomForest = false;
Params.frc6.SparseMethod = 'frc';
Params.frc6.Rotate = false;
Params.frc6.nmix = 6;

for c = 1:length(Classifiers)
    Params.(Classifiers{c}).ntrees = ntrees;
    Params.(Classifiers{c}).Stratified = Stratified;
    Params.(Classifiers{c}).NWorkers = NWorkers;
end

for i = 1:length(dims)
    
    d = dims(i);
    fprintf('dimension %d\n',d)
    
    for c = 1:length(Classifiers)
        
        cl = Classifiers{c};
        
        if Params.(cl).RandomForest
            if d <= 5
                Params.(cl).mtry{i} = 1:d;
            else
                Params.(cl).mtry{i} = ceil(d.^[0 1/4 1/2 3/4 1]);
            end
        else
            if d <= 5
                Params.(cl).mtry{i} = [1:d ceil(d.^[1.5 2])];
            elseif d > 5 && d <= 10
                Params.(cl).mtry{i} = ceil(d.^[0 1/4 1/2 3/4 1 1.5 2]);
            else
                Params.(cl).mtry{i} = [ceil(d.^[0 1/4 1/2 3/4 1]) 5*d 10*d];
            end
        end
        
        if isempty(Params.(cl).nmix) || Params.(cl).nmix <= d
            for j = 1:length(Params.(cl).mtry{i})
                Predictions = cell(n.test,n.trainSets);
                for trainSet = 1:n.trainSets
                    Forest = rpclassificationforest(Params.(cl).ntrees,...
                        X.train{i}(:,:,trainSet),Y.train{i}(:,trainSet),...
                        'RandomForest',Params.(cl).RandomForest,...
                        'sparsemethod',Params.(cl).SparseMethod,...
                        'nvartosample',Params.(cl).mtry{i}(j),...
                        'nmix',Params.(cl).nmix,...
                        'rotate',Params.(cl).Rotate,...
                        'Stratified',Stratified,...
                        'NWorkers',NWorkers);

                    Predictions(:,trainSet) = predict(Forest,X.test{i});
                end
                GE.(cl){i}(j) = misclassification_rate(Predictions,Y.test{i});
                V.(cl){i}(j) = classifier_variance(Predictions);
                B.(cl){i}(j) = classifier_bias(Predictions,TestPosteriors{i});
            end        
        else
            GE.(cl){i} = [];
            V.(cl){i} = [];
            B.(cl){i} = [];
        end
    end
end

save([rerfPath 'RandomerForest/Results/Sparse_parity_bias_variance.mat'],...
    'GE','V','B','Params')