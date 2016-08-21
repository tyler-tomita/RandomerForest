function I = gini_impurity(Y)

    n = length(Y);
    ClassLabels = unique(Y);
    phats = hist(Y,ClassLabels)/n;
    I = sum(phats.*(1-phats));
    
end