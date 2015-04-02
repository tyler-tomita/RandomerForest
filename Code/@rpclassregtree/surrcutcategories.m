function c=surrcutcategories(t,j)
%SURRCUTCATEGORIES Categories used for surrogate splits in decision tree.
%   C=SURRCUTCATEGORIES(T) returns an N-element cell array C of the
%   categories used for surrogate splits in the decision tree T, where N is
%   the number of nodes in the tree. For each node K, C{K} is a cell array.
%   The length of C{K} is equal to the number of surrogate predictors found
%   at this node. Every element of C{K} is either an empty string for a
%   continuous surrogate predictor or a two-element cell array with
%   categories for a categorical surrogate predictor. The first element of
%   this two-element cell array lists categories assigned to the left child
%   by this surrogate split and the second element of this two-element cell
%   array lists categories assigned to the right child by this surrogate
%   split. The order of the surrogate split variables at each node is
%   matched to the order of variables returned by SURRCUTVAR. The
%   optimal-split variable at this node is not included. For non-branch
%   (leaf) nodes, C contains an empty cell.
%
%   C=SURRCUTCATEGORIES(T,J) takes an array J of node numbers and returns
%   the categories for the specified nodes.
%
%   See also CLASSREGTREE, CLASSREGTREE/CUTCATEGORIES,
%   CLASSREGTREE/SURRCUTVAR, CLASSREGTREE/SURRCUTPOINT, CLASSREGTREE/SURRCUTTYPE. 

%   Copyright 2010 The MathWorks, Inc. 


if nargin>=2 && ~validatenodes(t,j)
    error(message('stats:classregtree:surrcutcategories:InvalidNode'));
end

if nargin<2
    j = 1:length(t.surrvar);
end

N = numel(j);
c = repmat({{}},N,1);

for i=1:N
    node = j(i);
    var = t.surrvar{node};
    if ~isempty(var)
        cut = t.surrcut{node};
        thisc = cell(1,numel(var));
        thisc(var>0) = {''};
        thisc(var<0) = cut(var<0);
        c{i} = thisc;
    end
end
