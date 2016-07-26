function Indices = stratified_cv(Y,K)
% STRATIFIED_CV generates cross-validation indices stratified by class
% label to ensure that each fold has a representative proportion of
% examples from each class. Indices is a vector of length n in which the
% ith element indicates which fold (1 through K) the ith training example
% belongs to.

Labels = unique(Y);
Indices = zeros(size(Y));

for i = 1:length(Labels)
    n = sum(Y==Labels(i));
    Indices(Y==Labels(i)) = crossvalind('Kfold',n,K);
end

end