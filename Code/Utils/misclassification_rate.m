function L = misclassification_rate(Predictions,Y,Weighted,Average)
% MISCLASSIFICATION_RATE Proportion of predictions made by a set of
% classifiers on a set of observations that do not match the true class
% label
%
%   L = MISCLASSIFICATION_RATE(PREDICTIONS,Y,WEIGHTED) returns the misclassification
%   rate of a set h1,...,hN classifiers on a set x1,...,xn of observations.
%   Rows of the n-by-N cell array PREDICTIONS correspond to observations
%   and columns correspond to classifiers. That is, PREDICTIONS(i,j) is the
%   predicted class label made by the jth classifier on the ith
%   observation. Y is an n-by-1 cell array of true class labels. L is a
%   scalar value if AVERAGE is true, else it's a N-vector. ISWEIGHT is a logical value. Specifying true takes into
%   account whether classifiers have different numbers of predictions,
%   which is often the case when the set of classifiers are trained with
%   bagging.

if ~exist('Weighted')
    Weighted = true;
else
    if ~islogical(Weighted)
        Weighted = true;
    end
end

if ~exist('Average')
    Average = true;
else
    if ~islogical(Average)
        Average = true;
    end
end

[~,N] = size(Predictions);
OOB = ~cellfun(@isempty,Predictions);
Wrong = NaN(size(Predictions));

parfor cl = 1:N
    oobidx = OOB(:,cl);
    Pred_cl = Predictions(:,cl);
    Wrong_cl = Wrong(:,cl);
    Wrong_cl(oobidx) = ~strcmp(Pred_cl(oobidx),Y(oobidx));
    Wrong(:,cl) = Wrong_cl;
%     Wrong(oobidx,cl) = ~strcmp(Predictions(OOB(:,cl),cl),Y(OOB(:,cl)));
end

if Average
    if Weighted
        L = nanmean(Wrong(:));
    else
        L = nanmean(nanmean(Wrong));
    end
else
    L = nanmean(Wrong);
end