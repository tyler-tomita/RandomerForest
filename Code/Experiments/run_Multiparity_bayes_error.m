close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Multiparity_data
load Random_matrix_adjustment_factor

Labels = unique(Ytest{1});
nLabels = length(Labels);

Posteriors = zeros(ntest,length(Labels));

%% Define joint distribution of x as a Gaussian mixture of eight components

% X is the data
% Y is the class label in {1,2,3}
% Z is the component label in {1,...,8}, where first two components are
% class 1, second two components are class 2, and last four are class 3

Mu = [0,0,-0.5*ones(1,p-2);...
    1,1,-0.5*ones(1,p-2);...
    1,0,-0.5*ones(1,p-2);...
    0,1,-0.5*ones(1,p-2);...
    0,0,0.5*ones(1,p-2);...
    1,1,0.5*ones(1,p-2);...
    1,0,0.5*ones(1,p-2);...
    0,1,0.5*ones(1,p-2)];

Sigma = repmat(diag([1/64*ones(1,2),ones(1,p-2)]),1,1,8);

%Prior probability of the eight Gaussian mixture components
Priors = [1/6*ones(4,1);1/12*ones(4,1)];

% list label by component
LabelByComp = cellstr(num2str([ones(2,1);2*ones(2,1);3*ones(4,1)]));

nComp = length(LabelByComp);

%% Make predictions with Bayes classifier

% evaluate joint density of (X,Z) at each point


for j = 1:length(Xtest)

    f_xz = zeros(ntest,nComp);

    for i = 1:nComp
        f_xz(:,i) = Priors(i)*mvnpdf(Xtest{j},Mu(i,:),Sigma(:,:,i));
    end

    % evaluate joint density of (X,Y) at each point
    f_xy = zeros(ntest,nLabels);

    for i = 1:length(Labels)
        f_xy(:,i) = sum(f_xz(:,strcmp(LabelByComp,Labels{i})),2);
    end

    f_x = sum(f_xy,2);

    Posteriors = zeros(ntest,nLabels);

    for i = 1:length(Labels)
        Posteriors(:,i) = f_xy(:,i)./f_x;
    end

    Predictions{j} = predict_class(Posteriors,Labels);
end
    
BayesError = misclassification_rate([Predictions{1};Predictions{2};...
    Predictions{3};Predictions{4}],[Ytest{1};Ytest{2};...
    Ytest{3};Ytest{4}],false);

save([rerfPath 'RandomerForest/Results/Multiparity_bayes_error.mat'],...
    'BayesError')