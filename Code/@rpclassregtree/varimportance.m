function imp = varimportance(Tree)
%VARIMPORTANCE Compute embedded estimates of variable importance.
%   IMP = VARIMPORTANCE(T) computes estimates of variable importance for
%   tree T by summing changes in the risk due to splits on every variable
%   and dividing the sum by the number of tree nodes. If the tree is grown
%   without surrogate splits, this sum is taken over best splits found at
%   each branch node. If the tree is grown with surrogate splits, this sum
%   is taken over all splits at each branch node including surrogate
%   splits. The returned vector IMP has one element for each input variable
%   in the data used to train this tree. At each node, the risk is
%   estimated as node impurity if impurity was used to split nodes and node
%   error otherwise. This risk is weighted by the node probability.
%   Variable importance associated with this split is computed as the
%   difference between the risk for the parent node and the total risk for
%   the two children. 
%
%   See also CLASSREGTREE, CLASSREGTREE/RISK, CLASSREGTREE/SURRCUTVAR.

%   Copyright 2008 The MathWorks, Inc. 


% If surrogate splits were computed during tree construction, Tree.varimp
% has been filled for the optimal and all surrogate splits (including those
% with negative measure of association) at every branch. If surrogate splits 
% were not computed, need to compute Tree.varimp now using optimal splits
% only.
if isempty(Tree.varimp)

    % Init
    imp = zeros(1,Tree.npred);
    
    % If impurity has been computed, use it to estimate variable
    % importance. Otherwise use classification error.
    N = Tree.numnodes;
    if ~isempty(Tree.impurity)
        crit = 'impurity';
    else
        crit = 'error';
    end
    r     = risk(Tree,1:N,'criterion',crit);
    
    % Parent risk minus unsplit risk
    rdiff = risk(Tree,1:N,'criterion',crit,'mode','diff');
    
    % Sum risk changes over all splits
    vars = abs(Tree.var);
    for parent=1:N
        kids = Tree.children(parent,:);
        if length(kids)==2 && kids(1)>0 && kids(2)>0
            ivar = vars(parent);
            imp(ivar) = imp(ivar) ...
                + rdiff(parent) - r(kids(1)) - r(kids(2));
        end
    end
    nbranch = sum(Tree.isbranch);
    if nbranch>0
        imp = imp/nbranch;
    end
    
else
    imp = Tree.varimp;
end

end
