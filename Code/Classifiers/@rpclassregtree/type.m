function ttype=type(t)
%TYPE Type of tree.
%   TTYPE=TYPE(T) returns the type of the tree T.  TTYPE is 'regression'
%   for regression trees and 'classification' for classification trees.
%
%   See also CLASSREGTREE.

%   Copyright 2006 The MathWorks, Inc. 


ttype = t.method;
