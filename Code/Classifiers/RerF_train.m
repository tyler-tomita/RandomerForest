function [Forest,Params] = RerF_train(Xtrain,Ytrain,Params)
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

if ~isfield(Params,'RandomForest')
    Params.RandomForest = false;
end

if ~isfield(Params,'Method')
    Params.Method = 'sparse-adjusted';
end

if ~isfield(Params,'Rotate')
    Params.Rotate = false;
end

if ~isfield(Params,'mdiff')
    Params.mdiff = 'off';
end

if ~isfield(Params,'d')
    if Params.RandomForest
        if p <= 5
            Params.d = 1:p;
        else
            Params.d = ceil(p.^[1/4 1/2 3/4 1]);
        end
    else
        if p <= 5
            Params.d = [1:p ceil(p.^[1.5 2])];
        elseif p > 5 && p <= 10
            Params.d = ceil(p.^[0 1/4 1/2 3/4 1 1.5 2]);
        else
            Params.d = [ceil(p.^[0 1/4 1/2 3/4 1]) 5*p 10*p];
        end
    end
end

if ~isfield(Params,'dprime')
    load Random_matrix_adjustment_factor
    for i = 1:length(Params.d)
        Params.dprime(i) = ceil(Params.d(i)^(1/interp1(ps,slope,p)));
    end
end

if ~isfield(Params,'nmix')
    Params.nmix = 2;
end

if ~isfield(Params,'Robust')
    Params.Robust = false;
end

if ~isfield(Params,'Stratified')
    Params.Stratified = true;
end

if ~isfield(Params,'NWorkers')
    Params.NWorkers = 2;
end

%train classifier for all values of Params.d

for i = 1:length(Params.d)
    Forest{i} = rpclassificationforest(Xtrain{j},Ytrain{j},...
        'nTrees',Params.nTrees,...
        'RandomForest',Params.RandomForest,...
        'sparsemethod',Params.Method,...
        'rotate',Params.Rotate,...
        'mdiff',Params.mdiff,...
        'nvartosample',Params.d(i),...
        'dprime',Params.dprime(i),...
        'nmix',Params.nmix,...
        'Robust',Params.Robust,...
        'Stratified',Params.Stratified,...
        'NWorkers',Params.NWorkers);
end