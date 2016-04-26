function err = oob_error(Predictions,Y,treenum)
    if nargin < 3
        treenum = 'last';
    end

    nrows = size(Y,1);
    ntrees = size(Predictions,2);
    predmat = NaN(nrows,ntrees);
    Labels = unique(Y);

    for j = 1:length(Labels)
        predmat(strcmp(Predictions,Labels{j})) = j;
    end

    if strcmp(treenum,'every')
        err = NaN(ntrees,1);
        parfor i = 1:ntrees
            ensemblePredictions = mode(predmat(:,1:i),2);
            missing = isnan(ensemblePredictions);
            ensemblePredictions = Labels(ensemblePredictions(~missing));
            wrong = ~strcmp(ensemblePredictions,Y(~missing));
            err(i) = mean(wrong);
        end
    else
        ensemblePredictions = mode(predmat,2);
        missing = isnan(ensemblePredictions);
        ensemblePredictions = Labels(ensemblePredictions(~missing));
        wrong = ~strcmp(ensemblePredictions,Y(~missing));
        err = mean(wrong);       
    end
end     %function oob_error