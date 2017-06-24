function [B,V,SE,VE,BE] = bias_variance(Predictions,Posteriors)
%BIAS_VARIANCE    Estimates the squared bias, variance, and Bayes error of a classifier given
%predictions for a set of classifiers each trained on a different training
%set. Bias and variance are defined according to equations (19-22) in
%"Variance and Bias for General Loss Functions" (James 2003)
%
%   [B,V,SE,VE] = BIAS_VARIANCE(PREDICTIONS,POSTERIORS) returns the estimated
%   squared bias and variance of a classifier given predictions on a set of observations
%   x1,...,xn for a set of classifiers h1,...,hN trained on training sets
%   D1,...,DN, respectively. PREDICTIONS is an n-by-N cell array of
%   predictions on n observations made by N classifiers. POSTERIORS is an
%   n-by-K matrix of posterior probabilities of being in each of K classes
%   given each of the n observations. That is, POSTERIORS(i,j) is
%   P(Y=yj|X=xi).

[n,N] = size(Predictions);
Labels = unique(Predictions);
Phats = zeros(n,length(Labels));
for i = 1:length(Labels)
    isY = false(size(Predictions));
    for cl = 1:N
        isY(:,cl) = strcmp(Predictions(:,cl),Labels{i});
    end
    Phats(:,i) = mean(isY,2);
end
[MaxPhat,PredictionIdx] = max(Phats,[],2);
Yhat = Labels(PredictionIdx);
[MaxPosterior,BayesIdx] = max(Posteriors,[],2);
Ybayes = Labels(BayesIdx);
B = mean(~strcmp(Yhat,Ybayes));
V = mean(1 - MaxPhat);
SE = mean(MaxPosterior - Posteriors(sub2ind(size(Posteriors),1:n,PredictionIdx)));
VE = mean(Posteriors(sub2ind(size(Posteriors),1:n,PredictionIdx)) - sum(Posteriors.*Phats,2));
BE = mean(1 - MaxPosterior);    % estimate of Bayes error
end