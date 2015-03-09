function xrank = map2rank(Xtrain,Xtrain_rank,xtest)
    xrank = zeros(1,size(Xtrain,2));
    Yrank = passtorank(cat(1,Xtrain,xtest));
    Yrank = Yrank(end,:);
    for i = 1:length(xtest)
        if Yrank(i) == 1
            xrank(i) = 1;
        elseif Yrank(i) > size(Xtrain,1)
            xrank(i) = size(Xtrain,1);
        else
            xrank(i) = (xtest(i) - Xtrain(Xtrain_rank(:,i)==Yrank(i)-1,i))/(Xtrain(Xtrain_rank(:,i)==Yrank(i),i) - Xtrain(Xtrain_rank(:,i)==Yrank(i)-1,i)) + Yrank(i) - 1;
        end
    end