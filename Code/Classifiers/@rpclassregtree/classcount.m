function p=classcount(t,j)
%CLASSCOUNT Class counts for tree nodes.
%   P=CLASSCOUNT(T) returns an N-by-M array P of class counts for the
%   nodes in the tree classification T, where N is the number of nodes
%   and M is the number of classes.  For any node number K, the class
%   counts P(K,:) are counts of observations (from the data used in
%   fitting the tree) from each class satisfying the conditions for node K.  
%
%   P=CLASSCOUNT(T,J) takes an array J of node numbers and returns the 
%   class counts for the specified nodes.
%
%   See also CLASSREGTREE, CLASSREGTREE/NUMNODES.

%   Copyright 2006-2007 The MathWorks, Inc. 


if isequal(t.method,'regression')
    error(message('stats:classregtree:classcount:NotClassification'));
end
if nargin>=2 && ~validatenodes(t,j)
    error(message('stats:classregtree:classcount:InvalidNode'));
end

if nargin<2
    p = t.classcount;
else
    p = t.classcount(j,:);
end
