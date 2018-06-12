function yhats = CCFPredict(xtrain,ytrain,xtest,nTrees,options)

    [CCF,preds,~,~] = genCCF(nTrees,xtrain,ytrain,false,options,xtest,false);
    yhats = CCF.classNames(preds);
end