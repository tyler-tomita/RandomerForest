function Xtest_rank = interpolate_rank(Xtrain,Xtest)
%Takes a test set and projects into the training data rank space
%For example suppose Xtrain is a 1-dimensional set of 3 points: [1;-2;5].
%Then passing to ranks gives Xtrain_rank = [2;1;3]. Now suppose you have a single point
%Xtest = 3. Since Xtest is half between Xtrain(1) and Xtrain(3), linearly
%interpolating gives its rank as being halfway between Xtrain_rank(1) and
%Xtrain_rank(3), or 2.5

%Xtrain: ntrain x d matrix of training data
%Xtest: ntest x d matrix of test data

    [ntest,d] = size(Xtest);
    ntrain = size(Xtrain,1);
    Xtest_rank = zeros(ntest,d);
    
    Xtrain_sort = sort(Xtrain);
    Xtrain_rank = tiedrank(Xtrain_sort);
    for i = 1:ntest
        for j = 1:d
            Xtest_idx = quickfind(Xtest(i,j),Xtrain_sort(:,j));
            if Xtest_idx == 0 || Xtest_idx == ntrain + 1;
                Xtest_rank(i,j) = Xtest_idx;
            else
                if rem(Xtest_idx,1) == 0
                    Xtest_rank(i,j) = Xtrain_rank(Xtest_idx,j);
                else
                    idx_below = floor(Xtest_idx);
                    idx_above = ceil(Xtest_idx);
                    X_below = Xtrain_rank(idx_below,j);
                    X_above = Xtrain_rank(idx_above,j);
                    Xtest_rank(i,j) = (Xtest_idx - idx_below)/(idx_above - idx_below)*(X_above - X_below) + X_below;
                end
            end
        end   
    end
end