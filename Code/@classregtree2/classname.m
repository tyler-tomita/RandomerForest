function cnames = classname(t,j)
%CLASSNAME Class names for classification decision tree.
%   CNAMES=CLASSNAME(T) returns a cell array of strings with class names
%   for this classification decision tree.
%
%   CNAMES=CLASSNAME(T,J) takes an array J of class numbers and returns the
%   class names for the specified numbers.
%
%   See also CLASSREGTREE.

%   Copyright 2010 The MathWorks, Inc. 


if nargin<2
    j = 1:length(t.classname);
end

cnames = t.classname(j);
end
