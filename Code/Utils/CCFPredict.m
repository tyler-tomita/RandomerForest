function yhats = CCFPredict(xtrain,ytrain,xtest,nTrees,options,iFeatureNum)

    if isempty(iFeatureNum)
        [CCF,preds,~,~] = genCCF(nTrees,xtrain,ytrain,false,options,xtest,false);
    else
        [CCF,preds,~,~] = genCCF(nTrees,xtrain,ytrain,false,options,xtest,false,iFeatureNum);
    end
    yhats = CCF.classNames(preds);
end