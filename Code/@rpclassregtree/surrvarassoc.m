function a=surrvarassoc(t,j)
%SURRVARASSOC Predictive measure of association for surrogate splits in decision tree.
%   A=SURRVARASSOC(T) returns an N-element cell array A of the predictive
%   measures of association for surrogate splits in the decision tree T,
%   where N is the number of nodes in the tree. For each node K, A{K} is a
%   numeric vector. The length of A{K} is equal to the number of surrogate
%   predictors found at this node. Every element of A{K} gives the
%   predictive measure of association between the optimal split and this
%   surrogate split. The order of the surrogate split variables at each
%   node is matched to the order of variables returned by SURRCUTVAR. The
%   optimal-split variable at this node is not included. For non-branch
%   (leaf) nodes, A contains an empty cell.
%
%   A=SURRVARASSOC(T,J) takes an array J of node numbers and returns the
%   predictive measure of association for the specified nodes.
%
%   See also CLASSREGTREE, CLASSREGTREE/SURRCUTVAR,
%   CLASSREGTREE/SURRCUTTYPE, CLASSREGTREE/SURRCUTCATEGORIES,
%   CLASSREGTREE/SURRCUTPOINT, CLASSREGTREE/SURRCUTFLIP.

%   Copyright 2010 The MathWorks, Inc.


if nargin>=2 && ~validatenodes(t,j)
    error(message('stats:classregtree:surrvarassoc:InvalidNode'));
end

if nargin<2
    j = 1:length(t.surrvar);
end

a = t.varassoc(j);
