function [importance,projections] = feature_importance(Forest,method,varargin)

    p = length(Forest.Tree{1}.rpm(:,1));
    Labels = Forest.classname;
    nSplitNodes = zeros(1,Forest.nTrees);
    
    for t = 1:Forest.nTrees
        nSplitNodes(t) = sum(Forest.Tree{t}.isbranch);
    end

    uniqueProjections = sparse(p,sum(nSplitNodes));

    idxstart = 1;
    for t = 1:Forest.nTrees
        idxend = idxstart + nSplitNodes(t) - 1;
        uniqueProjections(:,idxstart:idxend) = Forest.Tree{t}.rpm(:,Forest.Tree{t}.isbranch);
        idxstart = idxend + 1;
    end

    uniqueProjections = unique(uniqueProjections','rows')';
    nprojections = size(uniqueProjections,2);
    
    if strcmp(method,'gini')
        importance = total_gini(Forest,uniqueProjections');
        [importance,sortidx] = sort(importance,'descend');
        importance = (importance - importance(end))/(importance(1) - importance(end));
        projections = uniqueProjections(:,sortidx);
    elseif strcmp(method,'permutation')
        Xtrain = varargin{1};
        Ytrain = varargin{2};
        OOBError = varargin{3};
        n = length(Ytrain);
        importance = zeros(1,nprojections);
        for j = 1:nprojections
            fprintf('Computing importance of feature %d of %d\n',j,nprojections)
            projection = uniqueProjections(:,j);
            Ximp = full(Xtrain*projection);
            Ximp = Ximp(randperm(n));
            PermutedScores = rerf_oob_classprob(Forest,Xtrain,'last',true,projection,Ximp);
            PermutedPredictions = predict_class(PermutedScores,Labels);
            PermutedError = misclassification_rate(PermutedPredictions,Ytrain,false);
            importance(j) = PermutedError - OOBError;
        end
        [importance,sortidx] = sort(importance,'descend');
        importance = (importance - importance(end))/(importance(1) - importance(end));
        projections = uniqueProjections(:,sortidx);
    end
end

% function GiniImportance = total_gini(Forest,projection_imp)
%     GiniImportance = zeros(Forest.nTrees,1);
%     Trees = Forest.Tree;
%     parfor t = 1:Forest.nTrees
%         tree = Trees{t};
%         internalnodes = tree.node(tree.var ~= 0);
%         internalnodes = internalnodes';
%         leafnodes = tree.node(tree.var == 0);
%         leafnodes = leafnodes';
%         for node = internalnodes
%             projection = tree.rpm(:,node);
%             if isequal(projection,projection_imp)
%                 ch = tree.children(node);
%                 GiniImportance(t) = GiniImportance(t) + tree.impurity(node)*tree.nodesize(node) - ...
%                     (tree.impurity(ch(1))*tree.nodesize(ch(1)) + ...
%                     tree.impurity(ch(2))*tree.nodesize(ch(2)));
%             end
%         end
%     end
%     GiniImportance = sum(GiniImportance);
% end

function GiniImportance = total_gini(Forest,projections)
    % each row of projections corresponds to a particular projection
    nprojections = size(projections,1);
    GiniImportance = zeros(Forest.nTrees,nprojections);
    Trees = Forest.Tree;
    parfor t = 1:Forest.nTrees
        tree = Trees{t};
        internalnodes = tree.node(tree.var ~= 0);
        internalnodes = internalnodes';
        leafnodes = tree.node(tree.var == 0);
        leafnodes = leafnodes';
        for node = internalnodes
            projection = tree.rpm(:,node)';
            projectionIdx = find(ismember(projections,projection,'rows'));
            ch = tree.children(node);
            dI = tree.impurity(node)*tree.nodesize(node) - ...
                    (tree.impurity(ch(1))*tree.nodesize(ch(1)) + ...
                    tree.impurity(ch(2))*tree.nodesize(ch(2)));
            deltaImpurity = [zeros(1,projectionIdx-1),dI,zeros(1,nprojections - projectionIdx)];            
            GiniImportance(t,:) = GiniImportance(t,:) + deltaImpurity;
        end
    end
    GiniImportance = sum(GiniImportance);
end