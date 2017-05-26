function M = randmat(d,k,method,varargin)
    if strcmp(method,'binary')
        rho = varargin{1};
        M = zeros(d,k);
        nnzsPerCol = round(k*d*rho);
        nzs=randperm(d*k,nnzsPerCol);
        npos = rand(nnzsPerCol,1) > 0.5;
        M(nzs(npos))=1;
        M(nzs(~npos))=-1;
        isnz = M~=0;
        OnlyNz = isnz & repmat(sum(isnz)==1,d,1);
        M(OnlyNz) = 1;
        M = sparse(unique(M(:,any(M))','rows','stable')');
    elseif strcmp(method,'continuous')
        rho = varargin{1};
        M = zeros(d,k);
        nnzsPerCol = round(k*d*rho);
        nzs=randperm(d*k,nnzsPerCol);
        M(nzs) = rand(1,nnzsPerCol)*2 - 1;
        isnz = M~=0;
        OnlyNz = isnz & repmat(sum(isnz)==1,d,1);
        M(OnlyNz) = 1;
        M = sparse(unique(M(:,any(M))','rows','stable')');
    elseif strcmp(method,'binary-adjusted')
        rho = varargin{1};
        kk = varargin{3};
        M = zeros(d,kk);
        nnzsPerCol = round(kk*d*rho);
        nzs=randperm(d*kk,nnzsPerCol);
        npos = rand(nnzsPerCol,1) > 0.5;
        M(nzs(npos))=1;
        M(nzs(~npos))=-1;
        isnz = M~=0;
        OnlyNz = isnz & repmat(sum(isnz)==1,d,1);
        M(OnlyNz) = 1;
        M = unique(M(:,any(M))','rows','stable')';
        M = M(:,1:min(k,size(M,2)));
        M = sparse(M);
    elseif strcmp(method,'continuous-adjusted')
        rho = varargin{1};
        kk = varargin{3};
        M = zeros(d,kk);
        nnzsPerCol = round(kk*d*rho);
        nzs=randperm(d*kk,nnzsPerCol);
        M(nzs) = rand(1,nnzsPerCol)*2 - 1;
        isnz = M~=0;
        OnlyNz = isnz & repmat(sum(isnz)==1,d,1);
        M(OnlyNz) = 1;
        M = unique(M(:,any(M))','rows','stable')';
        M = M(:,1:min(k,size(M,2)));
        M = sparse(M);
    elseif strcmp(method,'frc')
        nmix = varargin{2};
        M = zeros(d,k);
%         p = 1;
%         for i = 1:nmix-1
%             p = p*(d-i)/d;
%         end
%         kk = round(4*k/p);
%         go = true;
%         while go
%             idx = randi(d,nmix,kk);
%             idx = idx(:,all(diff(sort(idx)),1));
%             go = size(idx,2) < k;
%         end
%         idx = idx(:,1:k);
        idx = randperms(d,nmix,k);
        idx = (ndgrid(1:k,1:nmix)'-1)*d + idx;
        M(idx) = rand(1,nmix*k)*2 - 1;
        M = sparse(M);
    elseif strcmp(method,'uniform-nnzs-binary')
        nmix = varargin{2};
        min_nmix = min(nmix);
        max_nmix = max(nmix);
        M = zeros(d,k);
        idx = randperms(d,max_nmix,k);
        idx = (ndgrid(1:k,1:max_nmix)'-1)*d + idx;
        nnzsPerCol = nmix(randi(length(nmix),1,k));
        for i = 1:length(nmix)
            idx(nmix(i)+1:end,nnzsPerCol==nmix(i)) = NaN;
        end
        idx(isnan(idx(:))) = [];
        nnzsTotal = length(idx(:));
        ispos = rand(nnzsTotal,1) > 0.5;
        M(idx(ispos)) = 1;
        M(idx(~ispos)) = -1;
        M = sparse(M);
    elseif strcmp(method,'uniform-nnzs-continuous')
        nmix = varargin{2};
        min_nmix = min(nmix);
        max_nmix = max(nmix);
        M = zeros(d,k);
        idx = randperms(d,max_nmix,k);
        idx = (ndgrid(1:k,1:max_nmix)'-1)*d + idx;
        nnzsPerCol = nmix(randi(length(nmix),1,k));
        for i = 1:length(nmix)
            idx(nmix(i)+1:end,nnzsPerCol==nmix(i)) = NaN;
        end
        idx(isnan(idx(:))) = [];
        M(idx(:)) = rand(1,length(idx(:)))*2 - 1;
        M = sparse(M);
    elseif strcmp(method,'poisson')
        lambda = varargin{2};
        M = zeros(d,k);
        go = true;
        while go
            nnzsPerCol = poissrnd(lambda,1,k);
            go = ~any(nnzsPerCol);
        end
        nnzsPerCol(nnzsPerCol > d) = d;
        nmix = unique(nnzsPerCol);
        max_nmix = nmix(end);
        idx = randperms(d,max_nmix,k);
        idx = repmat(0:k-1,max_nmix,1)*d + idx;
        for i = 1:length(nmix)
            idx(nmix(i)+1:end,nnzsPerCol==nmix(i)) = NaN;
        end
        one_nnz_idx = idx(:,nnzsPerCol==1);
        if isempty(one_nnz_idx)
            one_nnz_idx = [];
        else
            one_nnz_idx(isnan(one_nnz_idx)) = [];
        end
        idx(isnan(idx(:))) = [];
        nnzsTotal = length(idx(:));
        ispos = rand(nnzsTotal,1) > 0.5;
        M(idx(~ispos)) = -1;        
        M([idx(ispos),one_nnz_idx]) = 1;
%         M = sparse(M);
        M = sparse(unique(M(:,any(M))','rows','stable')');
    end
end