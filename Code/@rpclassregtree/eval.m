function [idname,nodes,id]=eval(Tree,X,subtrees)
%EVAL Compute fitted value for decision tree applied to data.
%   YFIT = EVAL(T,X) takes a classification or regression tree T and a
%   matrix X of predictor values, and produces a vector YFIT of predicted
%   response values. For a regression tree, YFIT(J) is the fitted response
%   value for a point having the predictor values X(J,:).  For a
%   classification tree, YFIT(J) is the class into which the tree would
%   assign the point with data X(J,:).
%
%   YFIT = EVAL(T,X,SUBTREES) takes an additional vector SUBTREES of
%   pruning levels, with 0 representing the full, unpruned tree.  T must
%   include a pruning sequence as created by the CLASSREGTREE constructor or
%   the PRUNE method. If SUBTREES has K elements and X has N rows, then the
%   output YFIT is an N-by-K matrix, with the Ith column containing the
%   fitted values produced by the SUBTREES(I) subtree.  SUBTREES must be
%   sorted in ascending order. (To compute fitted values for a tree that is
%   not part of the optimal pruning sequence, first use PRUNE to prune the
%   tree.)
%
%   [YFIT,NODE] = EVAL(...) also returns an array NODE of the same size
%   as YFIT containing the node number assigned to each row of X.  The
%   VIEW method can display the node numbers for any node you select.
%
%   [YFIT,NODE,CNUM] = EVAL(...) is valid only for classification trees.
%   It returns a vector CNUM containing the predicted class numbers.
%
%   NaN values and values absent from the training data for categorical
%   predictors in the X matrix are treated as missing.  If EVAL encounters
%   a missing value when it attempts to evaluate the split rule at a branch
%   node and if the tree is grown without surrogate splits, EVAL cannot
%   determine whether to proceed to the left or right child node. Instead,
%   it sets the corresponding fitted value equal to the fitted value
%   assigned to the branch node. If the tree is grown with surrogate
%   splits, EVAL evaluates the surrogate split rules at this branch node
%   until it finds a surrogate split rule for a predictor whose value is
%   not missing. If predictors for all surrogate split rules have missing
%   values, EVAL sets the corresponding fitted value equal to the fitted
%   value assigned to this branch node.
%
%   For a tree T, the syntax [...]=T(X) or [...]=T(X,SUBTREES) also invokes
%   the EVAL method.
%
%   Example: Find predicted classifications for Fisher's iris data.
%      load fisheriris;
%      t = classregtree(meas, species);  % create decision tree
%      sfit = eval(t,meas);              % find assigned class names
%      mean(strcmp(sfit,species))        % proportion correctly classified
%
%   See also CLASSREGTREE, CLASSREGTREE/PRUNE, CLASSREGTREE/VIEW, CLASSREGTREE/TEST.

%   Copyright 2006-2007 The MathWorks, Inc. 


if ~isfloat(X)
    error(message('stats:classregtree:eval:BadData'));
end

[nr,nc] = size(X);
if nc~=Tree.npred
   error(message('stats:classregtree:eval:BadInput', Tree.npred));
end

if nargin<3
   subtrees = 0;
elseif numel(subtrees)>length(subtrees)
   error(message('stats:classregtree:eval:TooManySubtrees'));
elseif any(diff(subtrees)<0)
   error(message('stats:classregtree:eval:UnsortedSubtrees'));
end

if ~isempty(Tree.prunelist)
   prunelist = Tree.prunelist;
elseif ~isequal(subtrees,0)
   Tree = prune(Tree);
else
   prunelist = Inf(size(Tree.node));
end

usedpreds = 1:Tree.npred;
treevar = Tree.var;
if isempty(Tree.surrvar)
    usedpreds = usedpreds(ismember(usedpreds,abs(treevar)));
    [~,treevar] = ismember(abs(treevar),usedpreds);
    treevar = treevar.*sign(Tree.var);
end

if isempty(usedpreds)
    nodes = ones(size(X,1),numel(subtrees));
else
    nodes = classregtreeEval(X(:,usedpreds)',Tree.node,treevar,Tree.cut,Tree.children,...
        prunelist,subtrees,Tree.surrvar,Tree.surrcut,Tree.surrflip);
    nodes = nodes';
end

id = Tree.class(nodes);
if nr==1
    id = id(:)';
end

if isequal(Tree.method,'classification')
   idname = Tree.classname(id);
   if nr==1
       idname = idname(:)';
   end
else
   idname = id;
end
end
