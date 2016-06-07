function M = randmat(d,k,method,varargin)
    if strcmp(method,'dense')
        s = varargin{1};
        M = vec2mat(randsample([0 1],d*k,true,[1-1/s 1/s]),k);
        M(M==1) = rand(sum(M(:)==1),1)*2 - 1;
    %elseif strcmp(method,'dgaussian')
    %    s = varargin{1};
    %    M = vec2mat(randsample([-1 0 1],d*k,true,[1/(2*s) 1-1/s 1/(2*s)]),k);
    %    M(M==1) = sqrt((s-1))*randn(sum(M(:)==1),1) + 1;
    %    M(M==-1) = sqrt((s-1))*randn(sum(M(:)==-1),1) - 1;
    elseif strcmp(method,'fast')
        M = sparse(d,k);
        M(randsample(d*k,k,false)) = randsample([-1 1],k,true,[0.5 0.5]);
        M = M(:,any(M));
    %elseif strcmp(method,'jovo')
    %    M = sparse(d,k);
    %    stop = false;
    %    while ~stop
    %        nnz=unique(round(rand(k+5,1)*(k*d-1))+1);
    %        nnzs = length(nnz);
    %        stop = nnzs <= numel(M) && nnzs >= 2;
    %    end
    %    M(nnz(1:round(nnzs/2)))=1;
    %    M(nnz(round(nnzs/2)+1:end))=-1;
    %    M = M(:,any(M));
    elseif strcmp(method,'sparse')
        s = varargin{1};
        kk = round(k/(1-1/exp(1)));
        M = sparse(d,kk);
        nnzs = round(kk*d*s);
        nzs=randperm(d*kk,nnzs);
        npos = rand(nnzs,1) > 0.5;
        M(nzs(npos))=1;
        M(nzs(~npos))=-1;
        M = M(:,any(M));
        M = M(:,1:min(k,size(M,2)));
    elseif strcmp(method,'sparse-uniform')
        s = varargin{1};
        kk = round(k/(1-1/exp(1)));
        M = sparse(d,kk);
        nnzs = round(kk*d*s);
        nzs = randperm(d*kk,nnzs);
        M(nzs) = rand(1,nnzs)*2 - 1;
        M = M(:,any(M));
        M = M(:,1:min(k,size(M,2)));
    elseif strcmp(method,'frc')
        nmix = varargin{2};
        M = sparse(d,k);
        p = 1;
        for i = 1:nmix-1
            p = p*(d-i)/d;
        end
        kk = round(4*k/p);
        go = true;
        while go
            idx = randi(d,nmix,kk);
            idx = idx(:,all(diff(sort(idx)),1));
            go = size(idx,2) < k;
        end
        idx = idx(:,1:k);
        idx = (ndgrid(1:k,1:nmix)'-1)*d + idx;
        M(idx) = rand(1,nmix*k)*2 - 1;
    elseif strcmp(method,'uniform')
        M = sparse(d,k);
        stop = false;
        while ~stop
            nnzs=unique(round(rand(k+5,1)*(k*d-1))+1);
            stop = length(nnzs) <= numel(M) && length(nnzs) >= 2;
        end
        M(nnzs)=rand(1,length(nnzs))*2 - 1;
        M = M(:,any(M));
    else
        M = sparse(d,k);
        R = poissrnd(1,1,k);
        %R(R==0) = 1;
        R(R>d) = d;
        rowidx = arrayfun(@(n) randsample(d,n,false),R,'UniformOutput',false);
        for j = 1:k
            M(rowidx{j},j) = 1;
        end
        %M(:,all(M==0,1)) = [];
        if strcmp(method,'gaussian')
            binsample = randsample([-1 1],sum(R),true,[0.5 0.5]);
            M(M==1) = binsample;
            M(M==1) = sqrt((d-1))*randn(sum(M(:)==1),1) + 1;
            M(M==-1) = sqrt((d-1))*randn(sum(M(:)==-1),1) - 1;
        else
            binsample = randsample([-1 1],sum(R),true,[0.5 0.5]);
            M(M==1) = binsample;
        end
    end
end