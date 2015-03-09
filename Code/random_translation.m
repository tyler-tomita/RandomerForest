%generates random translation matrix
function T = random_translation(n,d,sigma)
    T = repmat(sigma*randn(1,d),n,1);
end