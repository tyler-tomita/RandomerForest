function v=surrcutpoint(t,j)
%SURRCUTPOINT Cutpoints used for surrogate splits in decision tree.
%   V=SURRCUTPOINT(T) returns an N-element cell array V of the numeric
%   values used for surrogate splits in the decision tree T, where N is the
%   number of nodes in the tree. For each node K, V{K} is a numeric vector.
%   The length of V{K} is equal to the number of surrogate predictors found
%   at this node. Every element of V{K} is either NaN for a categorical
%   surrogate predictor or a numeric cut for a continuous surrogate
%   predictor. For every surrogate split with a numeric cut C based on a
%   continuous predictor variable Z, the left child is chosen if Z<C and
%   SURRCUTFLIP for this surrogate split is +1, or if Z>=C and SURRCUTFLIP
%   for this surrogate split is -1. Similarly, the right child is chosen if
%   Z>=C and SURRCUTFLIP for this surrogate split is +1, or if Z<C and
%   SURRCUTFLIP for this surrogate split is -1. The order of the surrogate
%   split variables at each node is matched to the order of variables
%   returned by SURRCUTVAR. The optimal-split variable at this node is not
%   included. For non-branch (leaf) nodes, V contains an empty cell.
%
%   V=SURRCUTPOINT(T,J) takes an array J of node numbers and returns the
%   cutpoints for the specified nodes.
%
%   See also CLASSREGTREE, CLASSREGTREE/CUTPOINT, CLASSREGTREE/SURRCUTVAR,
%   CLASSREGTREE/SURRCUTCATEGORIES, CLASSREGTREE/SURRCUTTYPE,
%   CLASSREGTREE/SURRCUTFLIP.

%   Copyright 2010 The MathWorks, Inc. 


if nargin>=2 && ~validatenodes(t,j)
    error(message('stats:classregtree:surrcutpoint:InvalidNode'));
end

if nargin<2
    j = 1:length(t.surrvar);
end

N = numel(j);
v = repmat({[]},N,1);

for i=1:N
    node = j(i);
    var = t.surrvar{node};
    if ~isempty(var)
        cut = t.surrcut{node};
        thisv = NaN(1,numel(var));
        thisv(var>0) = cell2mat(cut(var>0));
        v{i} = thisv;
    end
end
