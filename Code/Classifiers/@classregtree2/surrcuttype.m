function c=surrcuttype(t,j)
%SURRCUTTYPE Types of surrogate splits used at branches in decision tree.
%   C=SURRCUTTYPE(T) returns an N-element cell array C indicating types of
%   surrogate splits at each node in the tree T, where N is the number of
%   nodes in the tree. For each node K, C{K} is a cell array with the types
%   of the surrogate split variables at this node. The variables are sorted
%   by the predictive measure of association with the optimal predictor in
%   the descending order, and only variables with the positive predictive
%   measure are included. The order of the surrogate split variables at
%   each node is matched to the order of variables returned by SURRCUTVAR.
%   The optimal-split variable at this node is not included. For non-branch
%   (leaf) nodes, C contains an empty cell. A surrogate split type can be
%   either 'continuous' if the cut is defined in the form Z<V for a
%   variable Z and cutpoint V or 'categorical' if the cut is defined by
%   whether Z takes a value in a set of categories.
%
%   C=SURRCUTTYPE(T,J) takes an array J of node numbers and returns the cut
%   types for the specified nodes.
%
%   See also CLASSREGTREE, CLASSREGTREE/NUMNODES, CLASSREGTREE/CUTTYPE,
%   CLASSREGTREE/SURRCUTVAR. 

%   Copyright 2010 The MathWorks, Inc.


if nargin>=2 && ~validatenodes(t,j)
    error(message('stats:classregtree:surrcuttype:InvalidNode'));
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
        thisc = cell(1,numel(var));
        thisc(var>0) = {'continuous'};
        thisc(var<0) = {'categorical'};
        c{i} = thisc;
    end
end
