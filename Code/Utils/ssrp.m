function [H,W] = ssrp(X,Y,nvartosample,k,m,nlayers)

[H,W] = srp(X,Y,nvartosample,k);
if nlayers > 1
    for layer = 2:nlayers
        H = tanh(H);
        k = m*floor(sqrt(k));
        nvartosample = round(sqrt(size(H,2)));
        [H,W] = srp(H,Y,nvartosample,k);
    end
end