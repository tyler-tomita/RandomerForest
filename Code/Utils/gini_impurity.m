function I = gini_impurity(Y) %#codegen

    assert(size(Y,1) <= 1000);
    assert(isa(Y,'double'));
    if iscell(Y)
        Y = grp2idx(Y);
    end
    n = length(Y);
    ClassLabels = unique(Y);
    phats = hist(Y,ClassLabels)/n;
    I = sum(phats.*(1-phats));
end