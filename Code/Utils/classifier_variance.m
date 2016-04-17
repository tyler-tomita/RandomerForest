function v = classifier_variance(Predictions)
%CLASSIFIER_VARIANCE    Estimated variance of a set of classifiers
%
%   V = CLASSIFIER_VARIANCE(PREDICTIONS) returns the variance of
%   predictions made by a set of classifiers h1,...,hN on a set of points
%   x1,...,xn. Rows of the n-by-N cell array PREDICTIONS correspond to
%   observations and columns correspond to classifiers. That is,
%   PREDICTIONS(i,j) is the predicted class label made by the jth
%   classifier on the ith observation. The variance V is a scalar value.

[n,N] = size(Predictions);
OOB = ~cellfun(@isempty,Predictions);
Labels = unique(Predictions(OOB));
Phats = NaN(n,length(Labels));

for i = 1:length(Labels)
    isY = NaN(size(Predictions));
    for cl = 1:N
        isY(OOB(:,cl),cl) = strcmp(Predictions(OOB(:,cl),cl),Labels{i});
    end
    Phats(:,i) = nanmean(isY,2);
end

v = 0.5*(1 - nanmean(nansum(Phats.^2,2)));