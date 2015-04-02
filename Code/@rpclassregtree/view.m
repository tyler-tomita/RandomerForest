function outfig=view(Tree,varargin)
%VIEW View classification or regression tree graphically.
%   VIEW(T) displays the decision tree T as computed by the CLASSREGTREE
%   constructor in a figure window.  Each branch in the tree is labeled
%   with its decision rule, and each terminal node is labeled with the
%   predicted value for that node.  Click on any node to get more
%   information about it.  The information displayed is specified by
%   leftmost pop-up menu at the top of the figure.
%
%   VIEW(T,'PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%
%      'names'       A cell array of names for the predictor variables,
%                    in the order in which they appear in the X matrix
%                    from which the tree was created (see TREEFIT)
%      'prunelevel'  Initial pruning level to display
%
%   For each branch node, the left child node corresponds to the points
%   that satisfy the condition, and the right child node corresponds to
%   the points that do not satisfy the condition.
%
%   Example:  Create and graph classification tree for Fisher's iris data.
%             The names are abbreviations for the column contents (sepal
%             length, sepal width, petal length, and petal width).
%      load fisheriris;
%      t = classregtree(meas, species);
%      view(t,'names',{'SL' 'SW' 'PL' 'PW'});
%
%   See also CLASSREGTREE, CLASSREGTREE/EVAL, CLASSREGTREE/TEST, CLASSREGTREE/PRUNE.

%   Copyright 2006-2010 The MathWorks, Inc.


% Process inputs
nvars = max(abs(Tree.var));
varnames = names(Tree);
okargs =   {'names' 'prunelevel'};
defaults = {varnames   0};
[varnames,curlevel] = internal.stats.parseArgs(okargs,defaults,varargin{:});
if isempty(varnames)
   varnames = cell(nvars,1);
   for j=1:nvars
      varnames{j} = sprintf('x%d',j);
   end
elseif length(varnames)<nvars
   error(message('stats:classregtree:view:InvalidNames'));
end

impl = classreg.learning.impl.TreeImpl.makeFromClassregtree(Tree);
fig = viewGraph(impl,Tree.classname,Tree.class,varnames,curlevel,...
    '/stats/classregtree.view.html');

if nargout>0
    outfig = fig;
end
