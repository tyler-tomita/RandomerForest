function [classname,classid] = nodeclass(t,j)
%NODECLASS Class values in nodes of classification decision tree.
%   NAME=NODECLASS(T) returns an N-element cell array with the names of the
%   most probable classes in each node of the tree T, where N is the number
%   of nodes in the tree. Every element of this array is a string equal to
%   one of the class names returned by CLASSNAME(T). For regression trees,
%   NODECLASS returns an empty cell array.
%
%   NAME=NODECLASS(T,J) takes an array J of node numbers and returns the
%   class names for the specified nodes.
%
%   [NAME,ID]=NODECLASS(...) also returns a numeric array with the class
%   index for each node. The class index is determined by the order of
%   classes returned by CLASSNAME.
%
%   See also CLASSREGTREE, CLASSREGTREE/NUMNODES, CLASSREGTREE/CLASSNAME.

%   Copyright 2010 The MathWorks, Inc. 


if nargin>=2 && ~validatenodes(t,j)
    error(message('stats:classregtree:nodeclass:InvalidNode'));
end

if nargin<2
    j = 1:numnodes(t);
end

if strcmp(t.method,'classification')
    classid = t.class(j);
    classname = t.classname(classid);
else
    classid = [];
    classname = {};
end
end
