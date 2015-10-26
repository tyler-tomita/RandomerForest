function Lhat = RerF_crossval(X,Y,parms,Ind)

Yhats = cell(size(Y));

folds = unique(Ind);

for k = folds'
    
    Xtrain = X(Ind~=k,:);
    Ytrain = Y(Ind~=k);
    Xtest = X(Ind==k,:);
    
    Yhats(Ind==k) = RerF_classify(Xtest,Xtrain,Ytrain,parms);
end

Lhat = sum(~strcmp(Y,Yhats))/length(Y);