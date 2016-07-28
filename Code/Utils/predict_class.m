function Predictions = predict(Scores,Labels)
% PREDICT predicts class labels given Scores. Scores is an n-by-n_classes
% matrix of scores ranging from 0 to 1. Labels is an array of length 
% n_classes specifying the set of unique class labels. Values of Scores in 
% the jth column represent the strength of belief that an observation is 
% associated with the class label specified by Labels(j).

[~,PredictionIdx] = max(Scores,[],2);

Predictions = Labels(PredictionIdx);

end