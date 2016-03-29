function dv = diversity(predictions)

    if iscell(predictions)
        [I,G] = grp2idx(unique(predictions(~cellfun(@isempty,predictions))));
        predmat = NaN(size(predictions));
        for i = 1:length(G)
            predmat(strcmp(predictions,G(i))) = I(i);
        end
        predictions = predmat;
    end

    MajorityVote = mode(predictions,2);
    d = NaN(size(predictions));
    
    for tree = 1:size(predictions,2)
        Voted = ~isnan(predictions(:,tree));
        d(Voted,tree) = predictions(Voted,tree) ~= MajorityVote(Voted);
    end
    
    dv = nanmean(d(:));
end