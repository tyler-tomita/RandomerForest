function ma=meansurrvarassoc(t,j)
%MEANSURRVARASSOC Mean predictive measure of association for surrogate splits in decision tree.
%   MA=MEANSURRVARASSOC(T) returns a p-by-p matrix with predictive measures
%   of association for p predictors. Element MA(I,J) is the predictive
%   measure of association averaged over surrogate splits on predictor J
%   for which predictor I is the optimal split predictor. This average is
%   computed by summing positive values of the predictive measure of
%   association over optimal splits on predictor I and surrogate splits on
%   predictor J and dividing by the total number of optimal splits on
%   predictor I, including splits for which the predictive measure of
%   association between predictors I and J is negative.
%
%   MA=MEANSURRVARASSOC(T,N) takes an array N of node numbers and returns
%   the predictive measure of association averaged over the specified
%   nodes.
%
%   See also CLASSREGTREE, CLASSREGTREE/SURRCUTVAR,
%   CLASSREGTREE/SURRCUTTYPE, CLASSREGTREE/SURRCUTCATEGORIES,
%   CLASSREGTREE/SURRCUTPOINT, CLASSREGTREE/SURRCUTFLIP,
%   CLASSREGTREE/SURRVARASSOC.

%   Copyright 2010 The MathWorks, Inc.


if nargin>=2 && ~validatenodes(t,j)
    error(message('stats:classregtree:meansurrvarassoc:InvalidNode'));
end

if nargin<2
    j = 1:length(t.surrvar);
end

% Keep only branch nodes
isbr = isbranch(t,j);
j = j(isbr);

% Init
N = numel(j);
p = length(names(t));
ma = zeros(p);
nsplit = zeros(p,1);

% Get the association matrix and lists of best and surrogate predictors
a = t.varassoc(j);
[~,bestvar] = cutvar(t,j);
[~,surrvar] = surrcutvar(t,j);

% Loop over optimal splits. Increase the split count by 1 for every node.
for i=1:N
    n = bestvar(i);
    nsplit(n) = nsplit(n) + 1;
    m = surrvar{i};
    if ~isempty(m)
        ma(n,m) = ma(n,m) + a{i};
    end
end

% Divide cumulative association by the number of optimal splits on each
% predictor
gt0 = nsplit>0;
ma(gt0,:) = bsxfun(@rdivide,ma(gt0,:),nsplit(gt0));

% Make sure the diagonal elements are 1
ma(1:p+1:end) = 1;
