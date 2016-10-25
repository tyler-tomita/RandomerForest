function [Forest,Params,TrainTime] = RerF_train(Xtrain,Ytrain,Params)
% Params is a structure specifying algorithm parameters with fields:
%   nTrees: number of trees
%   RandomForest: logical true indicates do regular random forest
%   Method: method for contructing random matrix in oblique forests
%   Rotate: logical true indicates randomly rotate bootstrapped data prior 
%           to inducing trees
%   mdiff: 'node' indicates compute mean difference vector at each split
%          node and evaluate candidate projection onto this vector
%   d: number of variables to sample in RF or width of random matrix
%   dprime: adjusted value of d when Method is 'sparse-adjusted'
%   Robust: rank transform data prior to training
%   Stratified: logical true indicates stratify bootstraps by class
%   NWorkers: number of workers for inducing trees in parallel

p = size(Xtrain,2);

% set defaults if empty

if ~isfield(Params,'nTrees')
    Params.nTrees = 500;
end

if ~isfield(Params,'ForestMethod')
    Params.ForestMethod = 'sparse-binary';
end

if ~isfield(Params,'Rotate')
    Params.Rotate = false;
end

if ~isfield(Params,'mdiff')
    Params.mdiff = 'off';
end

if ~isfield(Params,'d')
    if strcmp(Params.ForestMethod,'rf')
        if p <= 5
            Params.d = 1:p;
        else
            Params.d = ceil(p.^[1/4 1/2 3/4 1]);
        end
    else
        if p <= 5
            Params.d = [1:p ceil(p.^[1.5 2])];
        elseif p > 5 && p <= 10
            Params.d = ceil(p.^[1/4 1/2 3/4 1 1.5 2]);
        else
            Params.d = [ceil(p.^[1/4 1/2 3/4 1]) 5*p 10*p];
        end
    end
end

if ~isfield(Params,'dprime')
    load Random_matrix_adjustment_factor
    for i = 1:length(Params.d)
        Params.dprime(i) = ceil(Params.d(i)^(1/interp1(dims,slope,p)));
    end
end

if ~isfield(Params,'nmix')
    Params.nmix = 2;
end

if ~isfield(Params,'Rescale')
    Params.Rescale = 'off';
end

if ~isfield(Params,'Stratified')
    Params.Stratified = true;
end

if ~isfield(Params,'NWorkers')
    Params.NWorkers = 2;
end

TrainTime = NaN(1,length(Params.d));

%train classifier for all values of Params.d

if strcmp(Params.ForestMethod,'frc')
    for i = 1:length(Params.nmix)
        for j = 1:length(Params.d)
            tic;
            if Params.nmix(i) == 1 && Params.d(j) > p
                Forest{(i-1)*length(Params.d)+j} = [];
                TrainTime((i-1)*length(Params.d)+j) = NaN;
            else
            Forest{(i-1)*length(Params.d)+j} = rpclassificationforest(Xtrain,Ytrain,...
                'nTrees',Params.nTrees,...
                'ForestMethod',Params.ForestMethod,...
                'rotate',Params.Rotate,...
                'mdiff',Params.mdiff,...
                'nvartosample',Params.d(j),...
                'dprime',Params.dprime(j),...
                'nmix',Params.nmix(i),...
                'Rescale',Params.Rescale,...
                'Stratified',Params.Stratified,...
                'NWorkers',Params.NWorkers);
                TrainTime((i-1)*length(Params.d)+j) = toc;
            end
        end
    end
else
    for i = 1:length(Params.d)
        tic;
        Forest{i} = rpclassificationforest(Xtrain,Ytrain,...
            'nTrees',Params.nTrees,...
            'ForestMethod',Params.ForestMethod,...
            'rotate',Params.Rotate,...
            'mdiff',Params.mdiff,...
            'nvartosample',Params.d(i),...
            'dprime',Params.dprime(i),...
            'nmix',Params.nmix,...
            'Rescale',Params.Rescale,...
            'Stratified',Params.Stratified,...
            'NWorkers',Params.NWorkers);
        TrainTime(i) = toc;
    end
end