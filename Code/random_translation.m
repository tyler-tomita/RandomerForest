%generates random translation matrix
function T = random_translation(n,d,min,max)
    pdf = makedist('Uniform','Lower',min,'Upper',max);
    T = repmat(random(pdf,1,d),n,1);
end