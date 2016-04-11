function p =  gmm_class_post(X,Mu,Sigma,Priors)
%GMM_CLASS_POST    Posterior probability of belonging to each mixture
%component in a multivariate gaussian mixture model
%   P = GMM_CLASS_POST(X,MU,SIGMA) returns the posterior probabilities of
%   class memberships in a gaussian mixture model specified by mean MU,
%   covariance SIGMA, and uniform class priors. Rows of the N-by-D matrix X
%   correspond to observations or points, and columns correspond to
%   variables or coordinates. Mu is a K-by-D vector, in which rows
%   designate the means for each of the K classes. SIGMA is a D-by-D-by-K
%   matrix, in which pages of SIGMA designate the covariances for each of
%   the K classes. P is an N-by-K vector.
%
%   P = GMM_CLASS_POST(X,MU,SIGMA,PRIORS) returns the posterior
%   probabilities of class membership in a gaussian mixture model with mean
%   MU, covariance SIGMA, and class priors specified by PRIORS. PRIORS is a
%   1-by-K vector, in which elements specify the prior probabilities of an
%   observation belonging to each of the K classes.

if nargin<3
    error('gmm_class_post:TooFewInputs','Too few inputs');
elseif ndims(X)~=2
    error('gmm_class_post:InvalidData','X must be n-by-d');
end

[n,d] = size(X);
if d<1
    error('gmm_class_post:TooFewDimensions','X must have >= 1 columns');
end

[k,d2] = size(Mu);
[d3,d4,k2] = size(Sigma);
if d2~=d
    error('gmm_class_post:MeanSizeMismatch','Mu and X must have same number of dimensions');
elseif d3~=d4
    error('gmm_class_post:BadCovariance','Covariance matrix must be square');
elseif d3~=d
    error('gmm_class_post:CovSizeMismatch','Sigma and X must have same number of dimensions');
elseif k~=k2
    error('gmm_class_post:NumClassMismatch','Number of mixture components in Mu and Sigma do not match');
end

if nargin<4 || isempty(Priors)
    Priors = 1/k*ones(1,k);
else
    k3 = length(Priors);
    if k3~=k
        error('gmm_class_post:NumClassMismatch','Number of mixture components in Priors and Mu do not match');
    end
end

f_xy_inv = zeros(n,k);
log_f = zeros(n,k); % log class conditional density of X given Y

for l = 1:k
    log_f(:,l) = log_mvnpdf(X,Mu(l,:),Sigma(:,:,l));
end

for l = 1:k
    for ll = 1:k
        f_xy_inv(:,l) = f_xy_inv(:,l) + exp(log_f(:,ll) + log(Priors(ll))...
            - log_f(:,l) - log(Priors(l)));
    end
    
end

p = 1./f_xy_inv;    % class posterior probabilities