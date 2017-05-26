function [importance,projections] = feature_importance(Forest,Xtrain,Ytrain,OOBError)
    [n,p] = size(Xtrain);
    Labels = Forest.classname;

    nSplitNodes = zeros(1,Forest.nTrees);
    for t = 1:Forest.nTrees
        nSplitNodes(t) = sum(Forest.Tree{t}.isbranch);
    end

    % p = length(Forest.Tree{1}.rpm(:,1));

    projections = sparse(p,sum(nSplitNodes));

    idxstart = 1;
    for t = 1:Forest.nTrees
        idxend = idxstart + nSplitNodes(t) - 1;
        projections(:,idxstart:idxend) = Forest.Tree{t}.rpm(:,Forest.Tree{t}.isbranch);
        idxstart = idxend + 1;
    end

    projections = unique(projections','rows')';
    nprojections = size(projections,2);
    
    importance = zeros(1,nprojections);

    for j = 1:nprojections
        fprintf('Computing importance of feature %d of %d\n',j,nprojections)
        projection = projections(:,j);
        Ximp = full(Xtrain*projection);
        Ximp = Ximp(randperm(n));
        PermutedScores = rerf_oob_classprob(Forest,Xtrain,'last',true,projection,Ximp);
        PermutedPredictions = predict_class(PermutedScores,Labels);
        PermutedError = misclassification_rate(PermutedPredictions,Ytrain,false);
        importance(j) = PermutedError - OOBError;
    end

    [importance,sortidx] = sort(importance,'descend');
    importance = (importance - importance(end))/(importance(1) - importance(end));
    projections = projections(:,sortidx);
end