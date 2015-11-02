function Yhats = RerF_classify(Xtest,Xtrain,Ytrain,parms)
%Trains discriminant analysis classifier and predicts class labels on test
%points

rerf = rpclassificationforest(parms.ntrees,Xtrain,Ytrain,'RandomForest',parms.RandomForest,...
    'sparsemethod',parms.sparsemethod,'nvartosample',parms.mtry,...
    'NWorkers',parms.NWorkers,'mdiff',parms.mdiff,'Stratified',parms.Stratified,...
    'Robust',parms.Robust,'splitcriterion',parms.splitcriterion);

Yhats = predict(rerf,Xtest);

end