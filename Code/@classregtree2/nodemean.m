function nm = nodemean(t,j)
%NODECLASS Mean values in nodes of regression decision tree.
%   NM=NODEMEAN(T) returns an N-element numeric array with mean values in
%   each node of the tree T, where N is the number of nodes in the tree.
%   Every element of this array is computed by averaging true Y values over
%   all observations in the node. For classification trees, NODEMEAN
%   returns an empty numeric array.
%
%   NM=NODEMEAN(T,J) takes an array J of node numbers and returns the
%   mean values for the specified nodes.
%
%   See also CLASSREGTREE, CLASSREGTREE/NUMNODES.

%   Copyright 2010 The MathWorks, Inc. 


if nargin>=2 && ~validatenodes(t,j)
    error(message('stats:classregtree:nodemean:InvalidNode'));
end

if nargin<2
    j = 1:numnodes(t);
end

if strcmp(t.method,'regression')
    nm = t.class(j);
else
    nm = [];
end
end
