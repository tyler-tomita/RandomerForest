function c=children(t,j)
%CHILDEREN Child nodes.
%   C=CHILDREN(T) returns an N-by-2 array C containing the numbers of the
%   child nodes for each node in the tree T, where N is the number of
%   nodes in the tree.  The children for leaf nodes (those with no children)
%   are returned as 0.
%
%   P=CHILDREN(T,J) takes an array J of node numbers and returns the children
%   for the specified nodes.
%
%   See also CLASSREGTREE, CLASSREGTREE/NUMNODES, CLASSREGTREE/PARENT.

%   Copyright 2006-2007 The MathWorks, Inc. 


if nargin>=2 && ~validatenodes(t,j)
    error(message('stats:classregtree:children:InvalidNode'));
end

if nargin<2
    c = t.children;
else
    c = t.children(j,:);
end

