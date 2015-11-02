function n=nodesize(t,j)
%NODESIZE Node size.
%   S=NODESIZE(T) returns an N-element vector S of the sizes of the
%   nodes in the tree T, where N is the number of nodes.  The size of a
%   node is defined as the number of observations from the data used to
%   create the tree that satisfy the conditions for the node.
%
%   S=NODESIZE(T,J) takes an array J of node numbers and returns the 
%   sizes for the specified nodes.
%
%   See also CLASSREGTREE, CLASSREGTREE/NUMNODES.

%   Copyright 2006 The MathWorks, Inc. 


if nargin>=2 && ~validatenodes(t,j)
    error(message('stats:classregtree:nodesize:InvalidNode'));
end

if nargin<2
    n = t.nodesize;
else
    n = t.nodesize(j,:);
end
