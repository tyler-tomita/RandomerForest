function [SplitValue,DeltaImpurity,Splits] = find_split_value(X,Y)

    if iscell(Y)
        Y = grp2idx(Y);
    end
    n = length(Y);
    [Xsorted,Idx] = sort(X);
    Ysorted = Y(Idx);
    isdiff = logical(diff(Ysorted));
    LowerIdx = [isdiff;false];
    UpperIdx = [false;isdiff];
    Splits = mean([Xsorted(LowerIdx),Xsorted(UpperIdx)],2);
    ParentImpurity = gini_impurity(Y);
    
    DeltaImpurity = zeros(length(Splits),1);
    
    for i = 1:length(Splits)
        yl = Ysorted(Xsorted<=Splits(i));
        yr = Ysorted(Xsorted>Splits(i));
        DeltaImpurity(i) = ParentImpurity - 1/n*(length(yl)*gini_impurity(yl) ...
            + length(yr)*gini_impurity(yr));
    end
    
    [~,BestIdx] = max(DeltaImpurity);
    if length(BestIdx) > 1
        BestIdx = randsample(BestIdx,1);
    end
    
    SplitValue = Splits(BestIdx);
    DeltaImpurity = DeltaImpurity(BestIdx);
    
end