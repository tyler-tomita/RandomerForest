function BestIdx = hp_optimize(CVError,AUC)
% HP_OPTIMIZE finds the index of the optimal hyperparameter value given
% either cross-validation errors, AUCs, or both

if ~isempty(CVError) && ~isempty(AUC)
    IsMinError = CVError==min(CVError);
    if sum(IsMinError)>1
        BestIdx = find(AUC==max(AUC(IsMinError)) & IsMinError);
    else
        BestIdx = find(IsMinError);
    end
elseif ~isempty(CVError) && isempty(AUC)
    BestIdx = find(CVError==min(CVError));
elseif isempty(CVError) && ~isempty(AUC)
    BestIdx = find(AUC==max(AUC));
end

end