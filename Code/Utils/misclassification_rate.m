function L = misclassification_rate(Predictions,Y)
% MISCLASSIFICATION_RATE Proportion of predictions made by a set of
% classifiers on a set of observations that do not match the true class
% label
%
%   L = MISCLASSIFICATION_RATE(PREDICTIONS,Y) returns the misclassification
%   rate of a set h1,...,hN classifiers on a set x1,...,xn of observations.
%   Rows of the n-by-N cell array PREDICTIONS correspond to observations
%   and columns correspond to classifiers. That is, PREDICTIONS(i,j) is the
%   predicted class label made by the jth classifier on the ith
%   observation. Y is an n-by-1 cell array of true class labels. L is a
%   scalar value.

[n,N] = size(Predictions);
OOB = ~cellfun(@isempty,Predictions);
Wrong = NaN(size(Predictions));

for cl = 1:N
    Wrong(OOB(:,cl),cl) = ~strcmp(Predictions(OOB(:,cl),cl),Y(OOB(:,cl)));
end

L = nanmean(Wrong(:));