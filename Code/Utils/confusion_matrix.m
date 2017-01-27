function C = confusion_matrix(Predictions,Y,Labels)
% MISCLASSIFICATION_RATE Computes the confusion matrix by comparing each
% element of PREDICTIONS with the corresponding true value in Y. C(i,j) is
% the number of observations known to be in group i but predicted to be in
% group j. LABELS is the unique set of possible values that Y can take.

n = length(Y);
nClasses = length(Labels);

C = zeros(nClasses,nClasses);

for i = 1:nClasses
    is_i = strcmp(Y,Labels{i});
    for j = 1:nClasses
        is_j = strcmp(Predictions,Labels{j});
        C(i,j) = sum(is_i & is_j);
    end
end