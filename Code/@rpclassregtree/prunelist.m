function p = prunelist(t,j)
%PRUNELIST Pruning levels for decision tree nodes.
%   P=PRUNELIST(T) returns an N-element numeric vector with the pruning
%   levels in each node of the tree T, where N is the number of nodes.
%   When you call PRUNE(T,'level',LEVEL), nodes with the pruning levels
%   below LEVEL are pruned, and nodes with the pruning levels greater or
%   equal to this level are not pruned.
%
%   P=PRUNELIST(T,J) takes an array J of node numbers and returns the
%   pruning levels for the specified nodes.
%
%   See also CLASSREGTREE, CLASSREGTREE/NUMNODES, CLASSREGTREE/PRUNE.

%   Copyright 2010 The MathWorks, Inc.


if nargin>=2 && ~validatenodes(t,j)
    error(message('stats:classregtree:prunelist:InvalidNode'));
end

if nargin<2
    j = 1:numnodes(t);
end

if isempty(t.prunelist)
    p = [];
else
    p = t.prunelist(j);
end
end
