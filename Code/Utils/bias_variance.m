function [B,V] = bias_variance(Predictions,Posteriors)
%BIAS_VARIANCE    Estimates the squared bias and variance of a classifier given
%predictions for a set of classifiers each trained on a different training
%set. Bias and variance are defined according to equations (19)
%and (20) in "Variance and Bias for General Loss Functions" (James 2003)
%
%   [B,V] = BIAS_VARIANCE(PREDICTIONS,POSTERIORS) returns the estimated
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
[~,PredictionIdx] = max(Phats,[],2);
Yhat = Labels(PredictionIdx);
[~,PredictionIdx] = max(Posteriors,[],2);
Ybayes = Labels(PredictionIdx);
B = mean(~strcmp(Yhat,Ybayes));
V = mean(1 - max(Phats,[],2));
end