function n=numnodes(t)
%NUMNODES Number of nodes in tree.
%   N=NUMNODES(T) returns the number of nodes N in the tree T.
%
%   See also CLASSREGTREE.

%   Copyright 2006 The MathWorks, Inc. 


n = numel(t.node);
