function p=nodeprob(t,j,mode)
%NODEPROB Node probability.
%   P=NODEPROB(T) returns an N-element vector P of the probabilities of the
%   nodes in the tree T, where N is the number of nodes.  The probability
%   of a node is computed as the proportion of observations from the
%   original data that satisfy the conditions for the node.  For a
%   classification tree, this proportion is adjusted for any prior
%   probabilities assigned to each class.
%
%   P=NODEPROB(T,J) takes an array J of node numbers and returns the 
%   probabilities for the specified nodes.
%
%   See also CLASSREGTREE, CLASSREGTREE/NUMNODES, CLASSREGTREE/NODESIZE.

%   Copyright 2006-2007 The MathWorks, Inc. 


if nargin>=2 && ~validatenodes(t,j)
    error(message('stats:classregtree:nodeprob:InvalidNode'));
end

if     nargin<2
    p = t.nodeprob;
    return;
end
j = j(:);

if     nargin<3
    p = t.nodeprob(j);
else
    if strcmpi(mode,'unsplit')
        p = zeros(numel(j),1);
        kids = t.children(j,:);
        isbr = all(kids>0,2);
        if any(isbr)
            p(isbr) = t.nodeprob(j(isbr)) - sum(t.nodeprob(kids(isbr,:)'),1)';
        end
    else
        error(message('stats:classregtree:nodeprob:BadMode'));
    end
end
end
