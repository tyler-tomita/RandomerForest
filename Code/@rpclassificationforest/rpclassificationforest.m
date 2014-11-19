classdef rpclassificationforest
    
    properties
        oobidx = {}; %indices of out of bag samples for each tree
        Tree = {};  %classification trees contained in the ensemble
        nTrees = []; %number of trees in the ensemble
        oobpredictions; %predicted class of oob samples for each tree...
        ooberror;   %fraction of training samples that were incorrectly classified
        classname;
    end
    
    methods
        function forest = rpclassificationforest(nTrees,X,Y,varargin)
            forest.oobpredictions = cell(size(X,1),nTrees);
            forest.ooberror = NaN(nTrees,1);
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
                        'SampleWithReplacement' 's'};
            defaults = {[]  []  'gdi'   []  []  1   sqrt(size(X,2))...
                        'off'    []  'off'    'classification'  1e-6    {}...
                        []  'off'   false  []  1    true   1};
            [Prior,Cost,Criterion,splitmin,minparent,minleaf,...
                nvartosample,Merge,categ,Prune,Method,qetoler,names,W,...
                surrogate,skipchecks,Stream,fboot,...
                SampleWithReplacement,s,~,extra] = ...
                internal.stats.parseArgs(okargs,defaults,varargin{:});
            nboot = ceil(fboot*length(Y));
            Tree = cell(nTrees,1);
            oobidx = cell(nTrees,1);
            parfor i = 1:nTrees
                sampleidx = 1:length(Y);
                ibidx = randsample(sampleidx,nboot,SampleWithReplacement);
                oobidx{i} = setdiff(sampleidx,ibidx);
                Tree{i} = rpclassregtree(X(ibidx,:),Y(ibidx,:),...
                    'priorprob',Prior,'cost',Cost,'splitcriterion',...
                    Criterion,'splitmin',splitmin,'minparent',...
                    minparent,'minleaf',minleaf,'nvartosample',...
                    nvartosample,'mergeleaves',Merge,'categorical',...
                    categ,'prune',Prune,'method',Method,'qetoler',...
                    qetoler,'names',names,'weights',W,'surrogate',...
                    surrogate,'skipchecks',skipchecks,'stream',Stream,...
                    's',s);
            end     %parallel loop over i
            forest.Tree = Tree;
            forest.oobidx = oobidx;
            forest.nTrees = length(forest.Tree);
        end     %function rpclassificationforest
        
        function forest = oobpredict(forest,X,Y)
            predmat = NaN(size(forest.oobpredictions));
            for i = 1:forest.nTrees
                forest.oobpredictions(forest.oobidx{i},i) = rptreepredict(forest.Tree{i},X(forest.oobidx{i},:));
            end
            for j = 1:length(forest.classname)
                predmat(strcmp(forest.oobpredictions,forest.classname{j})) = j;
            end
            for i = 1:forest.nTrees
                ensemblepredictions = mode(predmat(:,1:i),2);
                missing = isnan(ensemblepredictions);
                predictions = forest.classname(ensemblepredictions(~missing));
                wrong = ~strcmp(predictions,Y(~missing));
                forest.ooberror(i) = mean(wrong);
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