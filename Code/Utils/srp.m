function [H,W] = srp(X,Y,nvartosample,k)
%Semi-random projection

d = size(X,2);

W = zeros(d,k);

vars = logical(transpose(randerr(k,d,nvartosample)));
for j = 1:k
    Xr = X(:,vars(:,j));
    linclass = fitcdiscr(Xr,Y,'discrimType','pseudoLinear');
    W(vars(:,j),j) = linclass.Coeffs(1,2).Linear;
end

H = X*W;