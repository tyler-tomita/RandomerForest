function B = classifier_bias(Predictions,Posteriors)
%CLASSIFIER_BIAS    Estimates the squared bias of a classifier given
%predictions for a set of classifiers each trained on a different training
%set
%
%   B = CLASSIFIER_BIAS(PREDICTIONS,POSTERIORS) returns the estimated
%   squared bias of a classifier given predictions on a set of observations
%   x1,...,xn for a set of classifiers h1,...,hN trained on training sets
%   D1,...,DN, respectively. PREDICTIONS is an n-by-N cell array of
%   predictions on n observations made by N classifiers. POSTERIORS is an
%   n-by-K matrix of posterior probabilities of being in each of K classes
%   given each of the n observations. That is, POSTERIORS(i,j) is
%   P(Y=yj|X=xi).

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

B = 0.5*(nanmean(nansum((Posteriors - Phats).^2,2)));