classdef rpclassificationforest
    
    properties
        oobidx = {}; %indices of out of bag samples for each tree
        Tree = {};  %classification trees contained in the ensemble
        nTrees = []; %number of trees in the ensemble
        classname;
        RandomForest;
    end
    
    methods
        function forest = rpclassificationforest(nTrees,X,Y,varargin)
            if ~iscell(Y)
                Y = cellstr(num2str(Y));
            end
            forest.classname = unique(Y);
            forest = growTrees(forest,nTrees,X,Y,varargin{:});
        end     %class constructor
        
        function forest = growTrees(forest,nTrees,X,Y,varargin)
            okargs =   {'priorprob' 'cost'    'splitcriterion'  'splitmin'...
                        'minparent' 'minleaf'   'nvartosample'...
                        'mergeleaves'   'categorical' 'prune' 'method' ...
                        'qetoler'   'names'   'weights' 'surrogate'...
                        'skipchecks'    'stream'    'fboot'...
                        'SampleWithReplacement' 's' 'mdiff' 'sparsemethod' 'RandomForest'};
            defaults = {[]  []  'gdi'   []  []  1   ceil(size(X,2)^(2/3))...
                        'off'    []  'off'    'classification'  1e-6    {}...
                        []  'off'   false  []  1    true   1    'off'   'old' false};
            [Prior,Cost,Criterion,splitmin,minparent,minleaf,...
                nvartosample,Merge,categ,Prune,Method,qetoler,names,W,...
                surrogate,skipchecks,Stream,fboot,...
		  SampleWithReplacement,s,mdiff,sparsemethod,RandomForest,~,extra] = ...
                internal.stats.parseArgs(okargs,defaults,varargin{:});
            nboot = ceil(fboot*length(Y));
            Tree = cell(nTrees,1);
            oobidx = cell(nTrees,1);
            parfor i = 1:nTrees
                sampleidx = 1:length(Y);
                ibidx = randsample(sampleidx,nboot,SampleWithReplacement);
                oobidx{i} = setdiff(sampleidx,ibidx);
                if ~RandomForest
                    Tree{i} = rpclassregtree(X(ibidx,:),Y(ibidx,:),...
                        'priorprob',Prior,'cost',Cost,'splitcriterion',...
                        Criterion,'splitmin',splitmin,'minparent',...
                        minparent,'minleaf',minleaf,'nvartosample',...
                        nvartosample,'mergeleaves',Merge,'categorical',...
                        categ,'prune',Prune,'method',Method,'qetoler',...
                        qetoler,'names',names,'weights',W,'surrogate',...
                        surrogate,'skipchecks',skipchecks,'stream',Stream,...
                        's',s,'mdiff',mdiff,'sparsemethod',sparsemethod);
                else
                    Tree{i} = classregtree(X(ibidx,:),Y(ibidx,:),...
                        'priorprob',Prior,'cost',Cost,'splitcriterion',...
                        Criterion,'splitmin',splitmin,'minparent',...
                        minparent,'minleaf',minleaf,'nvartosample',...
                        nvartosample,'mergeleaves',Merge,'categorical',...
                        categ,'prune',Prune,'method',Method,'qetoler',...
                        qetoler,'names',names,'weights',W,'surrogate',...
                        surrogate,'skipchecks',skipchecks,'stream',Stream);
                end  
            end     %parallel loop over i
            forest.Tree = Tree;
            forest.oobidx = oobidx;
            forest.nTrees = length(forest.Tree);
            forest.RandomForest = RandomForest;
        end     %function rpclassificationforest
        
        function [err,varargout] = oobpredict(forest,X,Y,treenum)
            if nargin == 3
                treenum = 'last';
            end
            nrows = size(X,1);
            predmat = NaN(nrows,forest.nTrees);
            predcell = cell(nrows,forest.nTrees);
            OOBIndices = forest.oobidx;
            trees = forest.Tree;
            Labels = forest.classname;
            if ~forest.RandomForest
                for i = 1:forest.nTrees
                    pred_i = num2cell(NaN(nrows,1));
                    pred_i(OOBIndices{i}) = rptreepredict(trees{i},X(OOBIndices{i},:));
                    predcell(:,i) = pred_i;
                end
            else
                for i = 1:forest.nTrees
                    pred_i = num2cell(NaN(nrows,1));
                    pred_i(OOBIndices{i}) = eval(trees{i},X(OOBIndices{i},:));
                    predcell(:,i) = pred_i;
                end
            end
            for j = 1:length(forest.classname)
                predmat(strcmp(predcell,Labels{j})) = j;
            end
            if strcmp(treenum,'every')
                err = NaN(forest.nTrees,1);
                for i = 1:forest.nTrees
                    ensemblepredictions = mode(predmat(:,1:i),2);
                    missing = isnan(ensemblepredictions);
                    predictions = Labels(ensemblepredictions(~missing));
                    wrong = ~strcmp(predictions,Y(~missing));
                    err(i) = mean(wrong);
                end
            else
                ensemblepredictions = mode(predmat,2);
                missing = isnan(ensemblepredictions);
                predictions = Labels(ensemblepredictions(~missing));
                wrong = ~strcmp(predictions,Y(~missing));
                err = mean(wrong);         
            end
            if length(unique(Y)) == 2
                pos = num2str(max(str2num(char(Y))));
                neg = num2str(min(str2num(char(Y))));
                
                %varargout{1} = sum(strcmp(predictions(strcmp(pos,Y)),Y(~missing & strcmp(pos,Y))))/sum(strcmp(pos,Y));  %sensitivity
                %varargout{2} = sum(strcmp(predictions(strcmp(pos,Y)),Y(~missing & strcmp(pos,Y))))/sum(strcmp(pos,predictions));    %ppv
                %varargout{3} = sum(strcmp(predictions(strcmp(neg,Y)),Y(~missing & strcmp(neg,Y))))/sum(strcmp(neg,Y));  %specificity
                %varargout{4} = sum(strcmp(predictions(strcmp(neg,Y)),Y(~missing & strcmp(neg,Y))))/sum(strcmp(neg,predictions));    %npv
                varargout{1} = sum(strcmp(predictions(strcmp(pos,Y)),Y(~missing & strcmp(pos,Y)))); %tp
                varargout{2} = sum(~strcmp(predictions(strcmp(pos,Y)),Y(~missing & strcmp(pos,Y))));    %fn
                varargout{3} = sum(strcmp(predictions(strcmp(neg,Y)),Y(~missing & strcmp(neg,Y)))); %tn
                varargout{4} = sum(~strcmp(predictions(strcmp(neg,Y)),Y(~missing & strcmp(neg,Y))));    %fp
            end
        end     %function oobpredict
        
        function Y = predict(forest,X)
            predmat = NaN(size(X,1),forest.nTrees);
            YTree = cell(size(X,1),forest.nTrees);
            Tree = forest.Tree;
            parfor i = 1:forest.nTrees
                YTree(:,i) = rptreepredict(Tree{i},X);
            end
            classname = forest.classname;
            strcmp(YTree{1},classname{1});
            YTree{1};
            classname{1};
            for j = 1:length(classname)
                predmat(strcmp(YTree,classname{j})) = j;
            end
            predmat;
            ensemblepredictions = mode(predmat,2);
            missing = isnan(ensemblepredictions);
            Y = classname(ensemblepredictions(~missing));
        end     %function predict
    end     %methods
end     %classdef
