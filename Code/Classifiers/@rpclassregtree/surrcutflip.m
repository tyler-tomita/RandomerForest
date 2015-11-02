function v=surrcutflip(t,j)
%SURRCUTFLIP Numeric cutpoint assignments used for surrogate splits in decision tree.
%   V=SURRCUTFLIP(T) returns an N-element cell array V of the numeric cut
%   assignments used for surrogate splits in the decision tree T, where N
%   is the number of nodes in the tree. For each node K, V{K} is a numeric
%   vector. The length of V{K} is equal to the number of surrogate
%   predictors found at this node. Every element of V{K} is either zero for
%   a categorical surrogate predictor or a numeric cut assignment for a
%   continuous surrogate predictor. The numeric cut assignment can be
%   either -1 or +1. For every surrogate split with a numeric cut C based
%   on a continuous predictor variable Z, the left child is chosen if Z<C
%   and the cut assignment for this surrogate split is +1, or if Z>=C and
%   the cut assignment for this surrogate split is -1. Similarly, the right
%   child is chosen if Z>=C and the cut assignment for this surrogate split
%   is +1, or if Z<C and the cut assignment for this surrogate split is -1.
%   The order of the surrogate split variables at each node is matched to
%   the order of variables returned by SURRCUTVAR. The optimal-split
%   variable at this node is not included. For non-branch (leaf) nodes, V
%   contains an empty array.
%
%   V=SURRCUTFLIP(T,J) takes an array J of node numbers and returns the
%   cutpoint assignments for the specified nodes.
%
%   See also CLASSREGTREE, CLASSREGTREE/CUTPOINT, CLASSREGTREE/SURRCUTVAR,
%   CLASSREGTREE/SURRCUTCATEGORIES, CLASSREGTREE/SURRCUTTYPE,
%   CLASSREGTREE/SURRCUTPOINT. 

%   Copyright 2010 The MathWorks, Inc. 


if nargin>=2 && ~validatenodes(t,j)
    error(message('stats:classregtree:surrcutflip:InvalidNode'));
end

if nargin<2
    j = 1:length(t.surrvar);
end

v = t.surrflip(j);
