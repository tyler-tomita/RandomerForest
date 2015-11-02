function r=risk(t,j,varargin)
%RISK Node risk.
%   R=RISK(T) returns an N-element vector R of the risk of the nodes in the
%   tree T, where N is the number of nodes.  
%
%   R=RISK(T,J) takes an array J of node numbers and returns the risk
%   values for the specified nodes.
%
%   R=RISK(T,J,'criterion','error') returns risk vector R, where R(J) for
%   node J is the node error E(J) (classification error or mean squared
%   error for regression) weighted by the node probability P(J). 
%
%   R=RISK(T,J,'criterion','impurity') computes risk by using the node
%   impurity measure for each node instead of the node error. This option
%   is only valid for classification trees grown using impurity measures
%   such as Gini index or deviance. 
%
%   See also CLASSREGTREE, CLASSREGTREE/NODEERR, CLASSREGTREE/NODEPROB.

%   Copyright 2006-20097 The MathWorks, Inc. 


if nargin>=2 && ~validatenodes(t,j)
    error(message('stats:classregtree:risk:InvalidNode'));
end

if nargin<2
    j = 1:numnodes(t);
else
    j = j(:);
end

args = {'criterion'   'mode'};
defs = {'error'     'normal'};
[crit,mode] = internal.stats.parseArgs(args,defs,varargin{:});

if ~any(strncmpi(crit,{'error' 'impurity'},length(crit)))
    error(message('stats:classregtree:risk:BadCrit'));
end

if ~strncmpi(crit,'error',length(crit)) && isempty(t.impurity)
    error(message('stats:classregtree:risk:ImpurityNotAllowed'));
end

if ~any(strncmpi(mode,{'normal' 'unsplit' 'diff'},length(mode)))
    error(message('stats:classregtree:risk:BadMode'));
end

% Compute appropriate node probability
if     strncmpi(mode,'normal',length(mode))  % regular risk
    tprob = t.nodeprob(j);
elseif strncmpi(mode,'unsplit',length(mode)) % risk due to unsplit data
    tprob = nodeprob(t,j,'unsplit');
else                                         % difference: normal-unsplit
    tprob = t.nodeprob(j) - nodeprob(t,j,'unsplit');
end

% Compute risk
if strncmpi(crit,'error',length(crit))
    r = tprob.*t.nodeerr(j);
else
    r = tprob.*t.impurity(j);
end

end
