function [SplitValue,DeltaImpurity,Splits] = find_split_value(X,Y,Sorted,Yunique)
    if Sorted
        n = length(Y);
        Xunique = unique(X);
    %     isdiff = logical(diff(Y));
        isdiff = logical(diff(X)) & logical(diff(Y));
        LowerIdx = [isdiff;false];
        UpperIdx = [false;isdiff];
        Splits = mean([X(LowerIdx),X(UpperIdx)],2);
        RepeatIdx = ismember(Xunique,X(~logical(diff(X)) & logical(diff(Y))));
        LowerIdx = [RepeatIdx(2:end);false];
        UpperIdx = [false;RepeatIdx(2:end)];
        Splits = [Splits;mean([Xunique(LowerIdx),Xunique(UpperIdx)],2)];
        LowerIdx = [RepeatIdx(1:end-1);false];
        UpperIdx = [false;RepeatIdx(1:end-1)];
        Splits = [Splits;mean([Xunique(LowerIdx),Xunique(UpperIdx)],2)];
        Splits = unique(Splits);
        ParentImpurity = gini_impurity(Y,Yunique);

        DeltaImpurity = zeros(length(Splits),1);

        for i = 1:length(Splits)
            yl = Y(X<=Splits(i));   % this could be faster by using the SortIdx matrix
            yr = Y(X>Splits(i));
            DeltaImpurity(i) = ParentImpurity - 1/n*(length(yl)*gini_impurity(yl,Yunique) ...
                + length(yr)*gini_impurity(yr,Yunique));
        end

        [~,BestIdx] = max(DeltaImpurity);
        if length(BestIdx) > 1
            BestIdx = randsample(BestIdx,1);
        end

        SplitValue = Splits(BestIdx);
        DeltaImpurity = DeltaImpurity(BestIdx);
    else
        if iscell(Y)
            [Y,Yunique] = grp2idx(Y);
        end
        n = length(Y);
        [Xsorted,Idx] = sort(X);
        Xunique = unique(Xsorted);
        Ysorted = Y(Idx);
    %     isdiff = logical(diff(Ysorted));
        isdiff = logical(diff(Xsorted)) & logical(diff(Ysorted));
        LowerIdx = [isdiff;false];
        UpperIdx = [false;isdiff];
        Splits = mean([Xsorted(LowerIdx),Xsorted(UpperIdx)],2);
        RepeatIdx = ismember(Xunique,Xsorted(~logical(diff(Xsorted)) & logical(diff(Ysorted))));
        LowerIdx = [RepeatIdx(2:end);false];
        UpperIdx = [false;RepeatIdx(2:end)];
        Splits = [Splits;mean([Xunique(LowerIdx),Xunique(UpperIdx)],2)];
        LowerIdx = [RepeatIdx(1:end-1);false];
        UpperIdx = [false;RepeatIdx(1:end-1)];
        Splits = [Splits;mean([Xunique(LowerIdx),Xunique(UpperIdx)],2)];
        Splits = unique(Splits);
        ParentImpurity = gini_impurity(Y,Yunique);

        DeltaImpurity = zeros(length(Splits),1);

        for i = 1:length(Splits)
            yl = Ysorted(Xsorted<=Splits(i));
            yr = Ysorted(Xsorted>Splits(i));
            DeltaImpurity(i) = ParentImpurity - 1/n*(length(yl)*gini_impurity(yl,Yunique) ...
                + length(yr)*gini_impurity(yr,Yunique));
        end

        [~,BestIdx] = max(DeltaImpurity);
        if length(BestIdx) > 1
            BestIdx = randsample(BestIdx,1);
        end

        SplitValue = Splits(BestIdx);
        DeltaImpurity = DeltaImpurity(BestIdx);
    end
end