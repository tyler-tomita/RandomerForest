%generates random scaling matrix
function S = random_scaling(n,d,min,max)
    pdf = makedist('Uniform','Lower',min,'Upper',max);
    S = repmat(random(pdf,1,d),n,1);
end