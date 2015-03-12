function tf=isbranch(t,j)
%ISBRANCH Is node a branch node.
%   TF=ISBRANCH(T) returns an N-element logical vector TF that is true (1) for
%   each branch node and false (0) for each leaf node.
%
%   TF=ISBRANCH(T,J) takes an array J of node numbers and returns an array of
%   logical values for the specified nodes.
%
%   See also CLASSREGTREE, CLASSREGTREE/NUMNODES, CLASSREGTREE/CUTVAR.

%   Copyright 2006-2007 The MathWorks, Inc. 


if nargin>=2 && ~validatenodes(t,j)
    error(message('stats:classregtree:children:InvalidNode'));
end

if nargin<2
    tf = any(t.children~=0,2);
else
    tf = any(t.children(j,:)~=0,2);
end

