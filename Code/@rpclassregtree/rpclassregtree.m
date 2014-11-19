classdef rpclassregtree < classregtree

    properties(SetAccess='protected')
        rpm = {};
    end
    
    methods
        function a = rpclassregtree(X,Y,varargin)
            a = a@classregtree(X,Y);
            if nargin==1 && isa(X,'struct')
                a = struct2tree(a,X);            % convert from struct to tree
            else
                narginchk(2,Inf);
                a = rptreefit(a,X,Y,varargin{:});  % calls local version of treefit
                a.prunelist = zeros(0,1);
            end
        end % rpclassregtree constructor
        
        function Yhat = rptreepredict(Tree,X)
            Yhat = cell(size(X,1),1);
            noderows = cell(0,length(Tree.node));
            noderows{1} = 1:size(X,1);
            internalnodes = Tree.node(Tree.var ~= 0);
            internalnodes = internalnodes';
            leafnodes = Tree.node(Tree.var == 0);
            leafnodes = leafnodes';
            for node = internalnodes
                var = Tree.var(node);
                cut = Tree.cut{node};
                promat = Tree.rpm{node};
                Xpro = X(noderows{node},:)*promat(:,var);
                ch = Tree.children(node,:);
                noderows{ch(1)} = noderows{node}(Xpro < cut);
                noderows{ch(2)} = noderows{node}(Xpro >= cut);
            end
            for node = leafnodes
                Yhat(noderows{node}) = Tree.classname(Tree.class(node));
            end
        end     %function rptreepredict
    end     %methods block
end %classdef

function Tree=rptreefit(Tree,X,Y,varargin)

% Process inputs

X = sparse(X);  %reduce computational load for matrices with many zero elements

if isnumeric(Y)
   Method = 'regression';
else
   Method = 'classification';
end

okargs =   {'priorprob'   'cost'  'splitcriterion' ...
            'splitmin' 'minparent' 'minleaf' ...
            'nvartosample' 'mergeleaves' 'categorical' 'prune' 'method' ...
            'qetoler' 'names' 'weights' 'surrogate' 'skipchecks' ...
            'stream' 's'};
defaults = {[]            []      'gdi'                        ...
            []         []          1                          ...
            'all'          'on'          []            'off'    Method      ...
            1e-6      {}       []        'off'      false ...
            []  1};

[Prior,Cost,Criterion,splitmin,minparent,minleaf,...
    nvartosample,Merge,categ,Prune,Method,qetoler,names,W,surrogate,...
    skipchecks,Stream,s,~,extra] = ...
    internal.stats.parseArgs(okargs,defaults,varargin{:});

% For backwards compatibility. 'catidx' is a synonym for 'categorical'
for j=1:2:length(extra)
    if strncmp(extra{j},'catidx',length(extra{j}))
        categ = extra{j+1};
    else
        error(message('stats:classregtree:BadParamNameAndNotCatidx', extra{ j }));
    end
end

% Decode method
if ~ischar(Method) || isempty(Method) || ~(Method(1)=='c' || Method(1)=='r')
   error(message('stats:classregtree:BadMethod'));
elseif Method(1)=='c'
   Method = 'classification';
else
   Method = 'regression';
end

% Classification or regression?
doclass = strcmpi(Method(1),'c');

% Preprocess data
if isempty(W)
    W = ones(size(X,1),1);
end
if skipchecks
    if doclass
        % If checks are skipped, Y must be of type ClassLabel.
        [Y,cnames] = grp2idx(Y);
    end
else
    [X,Y,W,cnames] = classregtree.preparedata(X,Y,W,doclass);
end
[N,nvars] = size(X);

% Process variable names
if ~skipchecks
    names = classregtree.preparevars(names,nvars);
end

% Fill out criterion, class labels and matrix for classification
if doclass
   switch(Criterion)
    %                Criterion function   Is it an impurity measure?
    %                ------------------   --------------------------
    case 'gdi',      critfun = @gdi;      isimpurity = true;
    case 'twoing',   critfun = @twoing;   isimpurity = false;
    case 'deviance', critfun = @deviance; isimpurity = true;
    otherwise,     error(message('stats:classregtree:BadSplitCriterion'))
   end
   
   % Get binary matrix, C(i,j)==1 means point i is in class j
   nclasses = length(cnames);
   C = false(N,nclasses);
   C(sub2ind([N nclasses],(1:N)',Y)) = 1;   
   WC = bsxfun(@times,C,W);
   Wj = sum(WC,1);
else
   C = Y(:);
   isimpurity = [];
   critfun = [];
end

% Process prior and cost
if doclass
    % Get prior and cost
    if skipchecks
        if ~isstruct(Prior) || ~isstruct(Cost)
            error(message('stats:classregtree:BadPriorOrCostForSkipchecks'));
        end
        idx = getclassindex(cnames,Prior.group);
        Prior = Prior.prob(idx);
        idx = getclassindex(cnames,Cost.group);
        Cost = Cost.cost(idx,idx);
    else
        [Prior,Cost,removerows] = classregtree.priorandcost(Prior,Cost,cnames,Wj,Y);
        if any(removerows)
            X(removerows,:) = [];
            C(removerows,:) = [];
            WC(removerows,:) = [];
            Wj = sum(WC,1);
        end
        idx = Wj>0;
        W = sum(bsxfun(@times,WC(:,idx),Prior(idx)./Wj(idx)),2);
    end
    
    % Adjust prior to take costs into account
    % pratio is a multiplicative factor for class probabilities
    Cj = sum(Cost,2)';
    pratio = nclasses*Cj / sum(Cj);
    if ~isa(pratio,'double')
        pratio = double(pratio);
    end
else % regression
   pratio = 1;
end

% Clear WC to release memory
WC = [];

% Check and adjust node sizes.
if ~isempty(splitmin) && ~isempty(minparent)
    error(message('stats:classregtree:BothSplitminAndMinparentNotAllowed'));
end
if ~isempty(splitmin)
    if ~isnumeric(splitmin) || ~isscalar(splitmin)
        error(message('stats:classregtree:BadSplitmin'));
    end
    minparent = splitmin;
end
if ~isempty(minparent) && (~isnumeric(minparent) || ~isscalar(minparent))
    error(message('stats:classregtree:BadMinparent'));
end
if isempty(minparent)
    minparent = 10;
end
if ~isempty(minleaf) && (~isnumeric(minleaf) || ~isscalar(minleaf))
    error(message('stats:classregtree:BadMinleaf'));
end
if minleaf<1
    error(message('stats:classregtree:MinleafLessThanOne'));
end
if minparent<1
    error(message('stats:classregtree:MinparentLessThanOne'));
end
minparent = max(minparent,2*minleaf);

% Compute surrogate splits?
if ~strcmpi(surrogate,'on') && ~strcmpi(surrogate,'off')
    error(message('stats:classregtree:BadSurrogate'));
end
surrogate = strcmpi(surrogate,'on');

% Set the number of vars to be selected at random for splits
% Get number of features to sample
success = false;
if strcmpi(nvartosample,'all')
    nvartosample = nvars;
    success = true;
end
if isnumeric(nvartosample) && nvartosample>0 && nvartosample<=nvars
    nvartosample = ceil(nvartosample);
    success = true;
end
if ~success
    error(message('stats:classregtree:BadNumberOfRandomFeatures', nvars));
end
nusevars = nvartosample;

% Check prune parameter
if ~strcmpi(Prune,'on') && ~strcmpi(Prune,'off')
    error(message('stats:classregtree:BadPrune'));
end

% Check merge parameter
if ~strcmpi(Merge,'on') && ~strcmpi(Merge,'off')
    error(message('stats:classregtree:BadMergeLeaves'));
end

% Tree structure fields ([C] only for classification trees):
%  .method     method
%  .node       node number
%  .parent     parent node number
%  .class      class assignment for points in this node if treated as a leaf
%  .var        column j of X matrix to be split, or 0 for a leaf node,
%              or -j to treat column j as categorical
%  .cut        cutoff value for split (Xj<cutoff goes to left child node),
%              or a cell array of left and right values if var is negative
%  .children   matrix of child nodes (2 cols, 1st is left child)
%  .nodeprob   probability p(t) for this node
%  .nodeerr    resubstitution error estimate r(t) for this node
%  .nodesize   number of points at this node
%  .prunelist  list of indices that define pruned subtrees.  One entry per
%              node.  If prunelist(j)=k then, at the kth level of pruning,
%              the jth node becomes a leaf (or drops off the tree if its
%              parent also gets pruned).
%  .alpha      vector of complexity parameters for each pruning cut
%  .ntermnodes vector of terminal node counts for each pruning cut
%  .classprob  [C] vector of class probabilities
%  .classname  [C] names of each class
%  .classcount [C] count of members of each class
%  .nclasses   [C] number of classes
%  .cost       [C] misclassification cost

N = size(X,1);
Wtot = sum(W);
M = 2*ceil(N/minleaf)-1;% number of tree nodes for space reservation
nodenumber = zeros(M,1);
parent = zeros(M,1);
yfitnode = zeros(M,1);
cutvar = zeros(M,1);
cutpoint = cell(M,1);
children = zeros(M,2);
nodeprob = zeros(M,1);
resuberr = zeros(M,1);
nodesize = zeros(M,1);
rpm = cell(M,1);    %initialize cell array for storing proj matrices
if doclass
   classprob = zeros(M,nclasses);
   classcount = zeros(M,nclasses);
   if isimpurity
       impurity = zeros(M,1);
   end
   ybar = [];
end
iscat = false(nvars,1); iscat(categ) = 1;
nvarsplit = zeros(1,nvars);
varimp = [];
varassoc = [];
surrvar = {};
surrcut = {};
surrflip = {};
if surrogate
    varimp = zeros(1,nvars);
    varassoc = repmat({[]},M,1); % var associations
    surrvar = repmat({[]},M,1); % list of surrogate split vars for each node
    surrcut = repmat({{}},M,1); % list of surrogate split cuts for each node
    surrflip = repmat({[]},M,1);% list of flips: +1 do not flip, -1 flip
end

nodenumber(1) = 1;

assignednode = cell(M,1);% list of instances assigned to this node
assignednode{1} = 1:N;
nextunusednode = 2;

% Keep processing nodes until done
tnode = 1;
ybar = 0;
nclass = length(unique(Y));
%if nclass == 2
%    mu_diff = transpose(mean(X(Y==2,:)) - mean(X(Y==1,:)));
%end
while(tnode < nextunusednode)
   % Record information about this node
   noderows = assignednode{tnode};
   Nt = length(noderows);
   Cnode = C(noderows,:);
   Wnode = W(noderows);
   Wt = sum(Wnode);
   if doclass
      % Compute class probabilities and related statistics for this node
      Njt = sum(Cnode,1);    % number in class j at node t
      Pjandt = sum(bsxfun(@times,Cnode,Wnode),1);
      Pjgivent = Pjandt / sum(Pjandt);
      misclasscost = Pjgivent * Cost;
      [mincost,nodeclass] = min(misclasscost);
      yfitnode(tnode) = nodeclass;
      nodeprob(tnode) = Wt;
      classprob(tnode,:) = Pjgivent;
      classcount(tnode,:) = Njt;
      impure = sum(Pjgivent>0)>1;
      if isimpurity
          Pcorr = Pjgivent.*pratio;
          impurity(tnode) = feval(critfun,Pcorr/sum(Pcorr));
      end
   else
      % Compute variance and related statistics for this node
      ybar = sum(Cnode.*Wnode)/Wt;
      yfitnode(tnode) = ybar;
      nodeprob(tnode) = Wt/Wtot;
      sst = sum((Cnode-ybar).^2.*Wnode);% total sum of squares at this node
      mincost = sst / Wt;
      impure = (mincost > qetoler*resuberr(1));
   end
   bestcrit          = -Inf;
   nodesize(tnode)   = Nt;
   resuberr(tnode)   = mincost;
   cutvar(tnode)     = 0;
   cutpoint{tnode}   = 0;
   children(tnode,:) = 0;
   if surrogate
       varassoc(tnode) = {[]};
       surrvar(tnode) = {[]};
       surrcut(tnode) = {{}};
       surrflip(tnode) = {[]};
   end
   
   % Consider splitting this node
   if (Nt>=minparent) && impure      % split only large impure nodes
      Xnode = X(noderows,:);
      promat = srpmat(nvars,nusevars,s);    %random projection matrix
      %if nclass == 2
      %    promat = cat(2,mu_diff,promat);
      %end
      Xnode = Xnode*promat; %project Xnode onto random bases of promat
      bestvar = 0;
      bestcut = 0;
      
      % Find the best of all possible splits
      for jvar=1:size(Xnode,2)

         % Categorical variable?
         xcat = iscat(jvar);

         % Get rid of missing values and sort this variable
         idxnan = isnan(Xnode(:,jvar));
         idxnotnan = find(~idxnan);
         if isempty(idxnotnan)
             continue;
         end
         [x,idxsort] = sort(Xnode(idxnotnan,jvar));
         idx = idxnotnan(idxsort);
         c = Cnode(idx,:);
         w = Wnode(idx);
         
         % Downweight the impurity (for classification) or node mse (for
         % regression) by the fraction of observations that are being
         % split. Twoing already penalizes splits with low pL and pR.
         crit0U = 0;
         crit0  = 0;
         if doclass 
             if isimpurity % twoing crit does not need to be offset
                 % crit0U = P(t0-tU)*i(t0)
                 % crit0  = P(t0)*i(t0)
                 Pmis = sum(Wnode(idxnan));
                 crit0U = impurity(tnode)*(nodeprob(tnode)-Pmis);
                 crit0 = impurity(tnode)*nodeprob(tnode);
             end
         else
             % crit0U = P(t0-tU)*mse(t0)
             % crit0  = P(t0)*mse(t0)
             Pmis = sum(Wnode(idxnan));
             crit0U = resuberr(tnode)*(nodeprob(tnode)-Pmis);
             crit0 = resuberr(tnode)*nodeprob(tnode);
         end

         % Find optimal split for this variable
         [critval,cutval] = classregtreeRCcritval(full(x),doclass,c,w,pratio,...
             xcat,Criterion,bestcrit,double(crit0U),minleaf);
         
         % Change best split if this one is best so far
         if critval>bestcrit
            bestcrit = critval;
            bestvar = jvar;
            bestcut = cutval;
         end
      end

      % Split this node using the best rule found
      % Note: we have leftside==~rightside in the absence of NaN's
      if bestvar~=0
         nvarsplit(bestvar) = nvarsplit(bestvar)+1;
         x = Xnode(:,bestvar);
         
         % Send observations left or right
         if ~iscat(bestvar)
            cutvar(tnode) = bestvar;
            leftside = x<bestcut;
            rightside = x>=bestcut;
         else
            cutvar(tnode) = -bestvar;          % negative indicates cat. var. split
            leftside = ismember(x,bestcut{1});
            rightside = ismember(x,bestcut{2});
         end
         
         % Store split position, children, parent, and node number
         cutpoint{tnode} = bestcut;
         children(tnode,:) = nextunusednode + (0:1);
         nodenumber(nextunusednode+(0:1)) = nextunusednode+(0:1)';
         parent(nextunusednode+(0:1)) = tnode;
         rpm{tnode} = promat;
         
         %
         % Find surrogate splits
         %
         if surrogate
             % tsurrvar is an array with indices of
             %   surrogate vars (association with best var above zero)
             %   found for this split, excluding the best variable itself.
             %   These indices are the original var indices in input data.
             %
             % tvarassoc and tvarimp are numeric arrays with var
             %   associations (must be positive) and var importance values
             %   for these surrogate splits.
             %
             % tsurrcut is a cell array with surrogate split cuts, same
             %   convention as for cutpoint.
             %
             % tsurrflip is a numeric array with 0's for categorical
             %   surrogate splits and either -1 or +1 for numeric surrogate
             %   splits. -1 for a numeric splits means that left and right
             %   must be swapped, that is, leftside=x>=cut and
             %   rightside=x<cut for this surrogate split.
             %
             % tvarassoc, tsurrcut, and tsurflip have length
             %   numel(tsurrvar).
             %
             % tvarimp has length numel(varmap). These are variable
             % importance contributions from this branch node for *all*
             % predictors, not only surrogate predictors with positive
             % measure of association.
             %
             % tleftORright is a 2D array of size Nt-by-numel(tsurrvar)
             %   with surrogate split indices for observations: -1 if the
             %   surrogate split sends an observation left, +1 if it sends
             %   an observation right, and 0 if uncertain.
             [tvarassoc,tvarimp,tsurrvar,tsurrcut,tsurrflip,tleftORright] = ...
                 findsurrogate(Xnode,Cnode,Wnode,Wtot,doclass,isimpurity,critfun,...
                 varmap,iscat,bestvar,Cost,resuberr(tnode),pratio,crit0,...
                 leftside,rightside);
             
             % Update variable importance for cuts on best variable
             varimp(varmap) = varimp(varmap) + tvarimp;

             % Sort vars by their associations with the best var.
             [~,idxvarsort] = sort(tvarassoc,'descend');
             
             % Store surrogate cuts and if they need to be flipped
             surrcut(tnode) = {tsurrcut(idxvarsort)};
             surrflip(tnode) = {tsurrflip(idxvarsort)};
             
             % Store variables for surrogate splits.
             % For categorical vars, store negative indices.
             tsurrvar = tsurrvar(idxvarsort);
             tiscat = iscat(tsurrvar);
             tsurrvar(tiscat) = -tsurrvar(tiscat);
             surrvar(tnode) = {tsurrvar};
             
             % Store variable associations
             varassoc(tnode) = {tvarassoc(idxvarsort)};

             % Append lists of observations to be assigned to left and
             % right children
             for jmis=1:length(idxvarsort)
                 idxmissing = (1:Nt)';
                 idxmissing = idxmissing(~(leftside | rightside));
                 if isempty(idxmissing)
                     break;
                 else
                     surrmissing = tleftORright(idxmissing,idxvarsort(jmis));
                     leftside(idxmissing(surrmissing<0)) = true;
                     rightside(idxmissing(surrmissing>0)) = true;
                 end
             end             
         end
         
         % Assign observations for the next node
         assignednode{nextunusednode} = noderows(leftside);
         assignednode{nextunusednode+1} = noderows(rightside);
         
         % Update next node index
         nextunusednode = nextunusednode+2;
      end
   end
   tnode = tnode + 1;
end

topnode        = nextunusednode - 1;
Tree.method    = Method;
Tree.node      = nodenumber(1:topnode);
Tree.parent    = parent(1:topnode);
Tree.class     = yfitnode(1:topnode);
Tree.var       = cutvar(1:topnode);
Tree.cut       = cutpoint(1:topnode);
Tree.children  = children(1:topnode,:);
Tree.nodeprob  = nodeprob(1:topnode);
Tree.nodeerr   = resuberr(1:topnode);
Tree.nodesize  = nodesize(1:topnode);
Tree.npred     = nvars;
Tree.catcols   = categ;
Tree.names     = names;
Tree.minleaf   = minleaf;
Tree.minparent = minparent;
%if nclass == 2
%    Tree.nvartosample = nvartosample + 1;
%else
    Tree.nvartosample = nvartosample;
%end
Tree.mergeleaves = Merge;
Tree.nvarsplit = nvarsplit;
Tree.rpm = rpm(1:topnode);  %Store proj matrices in a structure field

if doclass
   Tree.prior     = Prior;
   Tree.nclasses  = nclasses;
   Tree.cost      = Cost;
   Tree.classprob = classprob(1:topnode,:);
   Tree.classcount= classcount(1:topnode,:);
   Tree.classname = cnames;
   if isimpurity
       Tree.impurity = impurity(1:topnode);
   end
   Tree.splitcriterion = Criterion;
else
    Tree.qetoler = qetoler;
end

% Get surrogate split info.
if ~isempty(surrvar)
    % Normalize var importance by the number of tree nodes.
    nbranch = sum(any(children(1:topnode,:)~=0,2));
    if nbranch>0
        Tree.varimp = varimp/nbranch;
    end

    % Throw away empty elements for surrogate splits.
    Tree.varassoc = varassoc(1:topnode);
    Tree.surrvar = surrvar(1:topnode);
    Tree.surrcut = surrcut(1:topnode);
    Tree.surrflip = surrflip(1:topnode);
end

if strcmpi(Merge,'on')
    Tree = mergeleaves(Tree); % merge leaves with same class
end

if strcmpi(Prune,'on')        % compute optimal pruning sequence if requested
   Tree = prune(Tree);
end
end

%----------------------------------------------------
function M = srpmat(d,k,s)
    M = rand(d,k);
    M(M <= 1/(2*s)) = -1*sqrt(s);
    M(M > 1/(2*s) & M <= 1-1/(2*s)) = 0;
    M(M > 1-1/(2*s)) = 1*sqrt(s);
    M = sparse(M);
end

%----------------------------------------------------
function v=gdi(p)
%GDI Gini diversity index

v=1-sum(p.^2,2);
end

%----------------------------------------------------
function v=twoing(Pleft, P1, Pright, P2)
%TWOING Twoing index

v = 0.25 * Pleft .* Pright .* sum(abs(P1-P2),2).^2;
end

%----------------------------------------------------
function v=deviance(p)
%DEVIANCE Deviance

v = -2 * sum(p .* log(max(p,eps(class(p)))), 2);
end

% --------------------------------------
function Tree = mergeleaves(Tree)
%MERGELEAVES Merge leaves that originate from the same parent node and give
% the sum of risk values greater or equal to the risk associated with the
% parent node.

N = length(Tree.node);
isleaf = (Tree.var==0)';   % no split variable implies leaf node
isntpruned = true(1,N);
doprune = false(1,N);
Risk = risk(Tree)';
unsplitRisk = risk(Tree,1:N,'mode','unsplit');
adjfactor = (1 - 100*eps(class(Risk)));

% Work up from the bottom of the tree
while(true)
   % Find ''twigs'' with two leaf children
   branches = find(~isleaf & isntpruned);
   twig = branches(sum(isleaf(Tree.children(branches,:)),2) == 2);
   if isempty(twig)
      break;            % must have just the root node left
   end
   
   % Find twigs to ''unsplit'' if the error of the twig is no larger
   % than the sum of the errors of the children
   Rtwig = Risk(twig);
   kids = Tree.children(twig,:);
   Rsplit = unsplitRisk(twig) + sum(Risk(kids),2);
   unsplit = Rsplit >= Rtwig'*adjfactor;
   if any(unsplit)
      % Mark children as pruned, and mark twig as now a leaf
      isntpruned(kids(unsplit,:)) = 0;
      twig = twig(unsplit);   % only these to be marked on next 2 lines
      isleaf(twig) = 1;
      doprune(twig) = 1;
   else
      break;
   end
end

% Remove splits that are useless
if any(doprune)
   Tree = prune(Tree,'nodes',find(doprune));
end
end

% ------------------------------------
function idx = getclassindex(cnames,g)
%GETCLASSINDEX Find indices for class names in another list of names
%   IDX = GETCLASSINDEX(CNAMES,G) takes a list CNAMES of class names
%   (such as the grouping variable values in the treefit or classify
%   function) and another list G of group names (as might be supplied
%   in the ''prior'' argument to those functions), and finds the indices
%   of the CNAMES names in the G list.  CNAMES should be a cell array
%   of strings.  G can be numbers, a string array, or a cell array of
%   strings

% Convert to common string form, whether input is char, cell, or numeric
if isnumeric(g)
   g = cellstr(strjust(num2str(g(:)), 'left'));
elseif ~iscell(g)
   g = cellstr(g);
end

% Look up each class in the grouping variable.
[~,idx] = ismember(cnames,g);
end

% ---------------------------------------
function [varassoc,varimp,surrvar,surrcut,surrflip,leftORright] = ...
    findsurrogate(Xnode,Cnode,Wnode,Wtot,doclass,isimpurity,critfun,...
    varmap,iscat,bestvar,Cost,tresuberr,pratio,crit0,leftside,rightside)
% Get number of vars and make default output
nvar = length(varmap);
N = size(Xnode,1);
varassoc = zeros(1,nvar);
varimp = zeros(1,nvar);
surrcut = cell(1,nvar);
surrvar = false(1,nvar);
surrflip = zeros(1,nvar);
leftORright = zeros(N,nvar);

% Total weight in this node
Wt = sum(Wnode);

% Left and right probabilities for the best split
pL = sum(Wnode(leftside))/Wt;
pR = sum(Wnode(rightside))/Wt;
minp = min(pL,pR);

% Loop over variables
for ivar=1:nvar
    % Get the predictor from the original data X
    jvar = varmap(ivar);
 
    % If best-split variable, assign left and right cases.
    % Best variable is not a surrogate variable but we need to compute
    % varimp for it.
    if jvar==bestvar
        leftORright(leftside,ivar)  = -1;
        leftORright(rightside,ivar) = +1;
    else
        %
        % Find the split that maximizes pLL+pRR
        %
        x = Xnode(:,jvar);        
        
        % If categorical variable, add every category to the side with
        % larger probability
        if iscat(jvar)
            [grp,~,grpval] = grp2idx(x);
            Ngrp = max(grp);
            if Ngrp<2
                continue;
            end
            sendgrpleft = false(Ngrp,1);
            for igrp=1:Ngrp
                tf = grp==igrp;
                Wleft = sum(Wnode(tf & leftside));
                Wright = sum(Wnode(tf & rightside));
                if Wleft<Wright
                    leftORright(tf,ivar) = +1;
                else
                    leftORright(tf,ivar) = -1;
                    sendgrpleft(igrp) = true;
                end
            end
            leftvals = grpval(sendgrpleft);
            rightvals = grpval(~sendgrpleft);
            if isempty(leftvals) || isempty(rightvals)
                continue;
            end
            leftvals = leftvals(:)';
            rightvals = rightvals(:)';
            pLL = sum(Wnode(leftORright(:,ivar)<0 & leftside))/Wt;
            pRR = sum(Wnode(leftORright(:,ivar)>0 & rightside))/Wt;
            if minp>1-pLL-pRR && pLL>0 && pRR>0
                surrvar(ivar) = true;
                surrcut{ivar} = {leftvals rightvals};
                varassoc(ivar) = (minp-(1-pLL-pRR)) / minp;
            end
            
        % If numeric variable, try all splits
        else
            % Find NaN's and sort
            idxnotnan = find(~isnan(x));
            if isempty(idxnotnan)
                continue;
            end
            [x,idxsorted] = sort(x(idxnotnan));
            idx = idxnotnan(idxsorted);
            
            % Determine if there's anything to split along this variable
            maxeps = max(eps(x(1)), eps(x(end)));
            if x(1)+maxeps > x(end)
                continue;
            end
            
            % Accept only splits on rows with distinct values
            idxdistinct = find(x(1:end-1) + ...
                max([eps(x(1:end-1)) eps(x(2:end))],[],2) < x(2:end));
            if isempty(idxdistinct)
                continue;
            end
            idxdistinct(end+1) = length(x);
            
            % Group into left and right using optimal split
            w = repmat(Wnode(idx)/Wt,1,2);
            w(rightside(idx),1) = 0;
            w(leftside(idx),2) = 0;
            w(~rightside(idx) & ~leftside(idx), :) = 0;
            w = cumsum(w,1);
            w = w(idxdistinct,:);
            x = x(idxdistinct);
            
            % Find split maximizing pLL+pRR
            [wLLandRRmax,i1] = ...
                max(w(1:end-1,1)+w(end,2)-w(1:end-1,2));
            [wLRandRLmax,i2] = ...
                max(w(end,1)-w(1:end-1,1)+w(1:end-1,2));
            if wLLandRRmax<wLRandRLmax
                surrflip(ivar) = -1;
                pLL = w(end,1)-w(i2,1);
                pRR = w(i2,2);
                cut = 0.5*(x(i2)+x(i2+1));
            else
                surrflip(ivar) = +1;
                pLL = w(i1,1);
                pRR = w(end,2)-w(i1,2);
                cut = 0.5*(x(i1)+x(i1+1));
            end
            x = Xnode(:,jvar);
            leftORright(x<cut,ivar)  = -surrflip(ivar);
            leftORright(x>=cut,ivar) = +surrflip(ivar);
            
            % Get association
            if minp>1-pLL-pRR && pLL>0 && pRR>0
                surrvar(ivar) = true;
                surrcut{ivar} = cut;
                varassoc(ivar) = (minp-(1-pLL-pRR)) / minp;
            end
        end
    end
        
    % Compute var importance
    sendleft = leftORright(:,ivar)<0;
    sendright = leftORright(:,ivar)>0;
    Cleft = Cnode(sendleft,:);
    Cright = Cnode(sendright,:);
    Wleft = Wnode(sendleft);
    Wright = Wnode(sendright);
    idxmiss = ~(sendleft | sendright);
    if doclass
        Pleft = sum(bsxfun(@times,Cleft,Wleft),1);
        Pright = sum(bsxfun(@times,Cright,Wright),1);
        if isimpurity
            Pleft = Pleft.*pratio;
            Pright = Pright.*pratio;
            varimp(ivar) = (1-sum(Wnode(idxmiss))/Wt)*crit0 ...
                - sum(Pleft)*feval(critfun,Pleft/sum(Pleft)) ...
                - sum(Pright)*feval(critfun,Pright/sum(Pright));
        else
            varimp(ivar) = Wt*tresuberr ...
                - min(Pleft*Cost) - min(Pright*Cost);
        end
    else
        ybarleft = sum(Wleft.*Cleft)/sum(Wleft);
        ybarright = sum(Wright.*Cright)/sum(Wright);
        varimp(ivar) = (1-sum(Wnode(idxmiss))/Wt)*crit0 ...
            - ( sum(Wleft.*(Cleft-ybarleft).^2) ...
            + sum(Wright.*(Cright-ybarright).^2) )/Wtot;
    end
end

% Return only values for surrogate split vars (satisfying varassoc>0).
% varimp is the only exception - it keeps values for all variables.
varassoc = varassoc(surrvar);
surrcut = surrcut(surrvar);
surrflip = surrflip(surrvar);
leftORright = leftORright(:,surrvar);
surrvar = varmap(surrvar);
end


% ------------------------------------
function treeobj=struct2tree(treeobj,S)
% Copy fields from structure S to tree object treeobj

% Look at all fields required for regression or classification trees
allfields = {'method'   'node'     'parent'   'class'   'var' ...
             'cut'      'children' 'nodeprob' 'nodeerr' ...
             'nodesize' 'npred'    'catcols'  ...
             'nclasses' 'prior'    'cost'     ...
             'classprob' 'classcount' 'classname'};
fn = fieldnames(S);
if ~ismember('method',fn) || ...
   (strcmpi(S.method,'classification') && ~all(ismember(allfields,fn))) || ...
   (strcmpi(S.method,'regression')     && ~all(ismember(allfields(1:12),fn)))
   error(message('stats:classregtree:BadTree'));
end
if strcmpi(S.method,'regression')
    nrequired = 12;
else
    nrequired = numel(allfields);
end
for j=1:nrequired
    fname = allfields{j};
    treeobj.(fname) = S.(fname);
end

% Look at optional fields
optionalfields = {'names' 'prunelist' 'alpha' 'ntermnodes' ...
    'impurity' 'prunecriterion' ...
    'minparent' 'minleaf' 'nvartosample' 'mergeleaves' ...
    'splitcriterion' 'qetoler' 'varassoc' 'varimp' 'nvarsplit' ...
    'surrvar' 'surrcut' 'surrflip' 'catsplit'};
for j=1:numel(optionalfields)
    fname = optionalfields{j};
    if isfield(S,fname)
        treeobj.(fname) = S.(fname);
    end
end
end
