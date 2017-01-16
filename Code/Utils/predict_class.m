function Predictions = predict(Scores,Labels)
% PREDICT predicts class labels given Scores. Scores is an n-by-nClasses
% matrix of scores ranging from 0 to 1. Labels is an array of length 
% nClasses specifying the set of unique class labels. Values of Scores in 
% the jth column represent the strength of belief that an observation is 
% associated with the class label specified by Labels(j).

n = size(Scores,1);

nanIdx = isnan(sum(Scores,2));

[~,PredictionIdx] = max(Scores(~nanIdx,:),[],2);

Predictions = cell(n,1);

Predictions(~nanIdx) = Labels(PredictionIdx);

end