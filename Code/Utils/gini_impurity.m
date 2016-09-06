function I = gini_impurity(Y)

    if iscell(Y)
        Y = grp2idx(Y);
    end
    n = length(Y);
    ClassLabels = unique(Y);
    phats = hist(Y,ClassLabels)/n;
    I = sum(phats.*(1-phats));
    
end