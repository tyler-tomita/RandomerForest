function P = randperms(N,K,M)
    [~,P] = sort(rand(N,M));
    P = P(1:K,:);
end