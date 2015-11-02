function display(a,objectname)
%@DATASET/DISPLAY Display a CLASSREGTREE object.
%   DISPLAY(T) prints the CLASSREGTREE object T.  DISPLAY is called when a
%   semicolon is not used to terminate a statement.
%
%   See also CLASSREGTREE, CLASSREGTREE/EVAL, CLASSREGTREE/TEST, CLASSREGTREE/PRUNE.

%   Copyright 2006-2007 The MathWorks, Inc. 


isLoose = strcmp(get(0,'FormatSpacing'),'loose');

if nargin<2
    objectname = inputname(1);
end
if isempty(objectname)
   objectname = 'ans';
end

if (isLoose), fprintf('\n'); end
fprintf('%s = \n', objectname);
disp(a)
