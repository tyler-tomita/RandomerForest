function I = gini_impurity(Y,Labels) %#codegen
    if iscell(Y)
        Y = grp2idx(Y);
    end
    n = length(Y);
    ClassCounts = zeros(1,length(Labels));
    for i = 1:length(Labels)-1
        ClassCounts(i) = sum(Y==Labels(i));
    end
    ClassCounts(end) = n - sum(ClassCounts(1:end-1));
    ClassProb = ClassCounts/n;
    I = sum(ClassProb.*(1-ClassProb));
end