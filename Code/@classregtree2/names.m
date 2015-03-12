function n = names(t,j)
%NAMES Variables used in decision tree.
%   N=NAMES(T) returns a cell array of strings with P elements for P
%   variables used to construct tree T.
%
%   N=NAMES(T,J) takes an array of indices J and returns the variables for
%   the specified indices.
%
%   See also CLASSREGTREE.

%   Copyright 2011 The MathWorks, Inc.


if isnumeric(t.names)
    n = strcat({'x'},num2str((1:t.names)','%-d'))';
else
    n = t.names;
end
if nargin>1
    n = n(j);
end
end
