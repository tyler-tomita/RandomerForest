%projection forest version 2

classdef rpclassificationforest2
    
    properties
        oobidx = {}; %indices of out of bag samples for each tree
        Trees = {};  %classification trees contained in the ensemble
        nTrees = []; %number of trees in the ensemble
        rpm = {};   %random projection matrices
        classname;
    end
    
    methods
        function forest = rpclassificationforest2(nTrees,X,Y,varargin)
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
                        'SampleWithReplacement' 's' 'mdiff' 'sparsemethod'...
                        'RandomForest'  'VarImp'};
            defaults = {[]  []  'gdi'   []  []  1   ceil(size(X,2)^(2/3))...
                        'off'    []  'off'    'classification'  1e-6    {}...
                        []  'off'   false  []  1    true   1    'off'   'old'...
                        false []};
            [Prior,Cost,Criterion,splitmin,minparent,minleaf,...
                nvartosample,Merge,categ,Prune,Method,qetoler,names,W,...
                surrogate,skipchecks,Stream,fboot,...
                SampleWithReplacement,s,mdiff,sparsemethod,RandomForest,VarImp,~,extra] = ...
                internal.stats.parseArgs(okargs,defaults,varargin{:});
            if isnumeric(nvartosample)
                nvartosample = ceil(nvartosample);
            elseif strcmp(nvartosample,'all')
                nvartosample = size(X,2);
            end
            nvars = size(X,2);
            nboot = ceil(fboot*length(Y));
            Tree = cell(nTrees,1);
            oobidx = cell(nTrees,1);
            if ~RandomForest
                promat = cell(nTrees,1);
                if strcmp(mdiff,'on')
                    Labels_str = unique(Y);
                    Labels = str2num(char(Labels_str));
                    K = length(Labels);
                    pairs = zeros(K-1,2);
                    npairs = K-1;
                    pairs(:,1) = 1:npairs;
                    pairs(:,2) = pairs(:,1) + 1;
                    %mu_diff = zeros(nvars,npairs);
                    %for i = 1:npairs
                    %    mu_diff(:,i) = transpose(mean(X(strcmp(Y,Labels_str(pairs(i,2))),:)) - mean(X(strcmp(Y,Labels_str(pairs(i,1))),:)));
                    %end
                    if isempty(VarImp)
                        Imp = ones(nvars,npairs);
                    else
                        b = TreeBagger(nTrees,X,Y,'NVarToSample',nvartosample,'OOBPred','on','OOBVarImp','on');
                        Imp = b.OOBPermutedVarDeltaError';
                        clear b
                        Imp(Imp<VarImp) = 0;
                        Imp = repmat(Imp,1,npairs);
                    end
                    parfor i = 1:nTrees
                        sampleidx = 1:length(Y);
                        ibidx = randsample(sampleidx,nboot,SampleWithReplacement);
                        oobidx{i} = setdiff(sampleidx,ibidx);
                        Xib = X(ibidx,:);
                        Yib = Y(ibidx);
                        mu_diff = zeros(nvars,npairs);
                        for j = 1:npairs
                            mu_diff(:,j) = transpose(mean(Xib(strcmp(Yib,Labels_str(pairs(j,2))),:)) - mean(Xib(strcmp(Yib,Labels_str(pairs(j,1))),:)));
                        end
                        mu_diff = mu_diff.*Imp;
                        promat{i} = srpmat(nvars,nvars,sparsemethod,s);    %random projection matrix
                        promat{i} = cat(2,mu_diff,promat{i});
                        Xpro = Xib*promat{i};
                        Tree{i} = classregtree2(Xpro,Yib,...
                            'priorprob',Prior,'cost',Cost,'splitcriterion',...
                            Criterion,'splitmin',splitmin,'minparent',...
                            minparent,'minleaf',minleaf,'nvartosample',...
                            nvartosample,'mergeleaves',Merge,'categorical',...
                            categ,'prune',Prune,'method',Method,'qetoler',...
                            qetoler,'names',names,'weights',W,'surrogate',...
                            surrogate,'skipchecks',skipchecks,'stream',Stream);
                    end     %parallel loop over i
                else
                    parfor i = 1:nTrees
                    sampleidx = 1:length(Y);
                    ibidx = randsample(sampleidx,nboot,SampleWithReplacement);
                    oobidx{i} = setdiff(sampleidx,ibidx);
                    promat{i} = srpmat(nvars,nvars,sparsemethod,s);    %random projection matrix
                    Xpro = X(ibidx,:)*promat{i};
                    Tree{i} = classregtree2(Xpro,Y(ibidx),...
                        'priorprob',Prior,'cost',Cost,'splitcriterion',...
                        Criterion,'splitmin',splitmin,'minparent',...
                        minparent,'minleaf',minleaf,'nvartosample',...
                        nvartosample,'mergeleaves',Merge,'categorical',...
                        categ,'prune',Prune,'method',Method,'qetoler',...
                        qetoler,'names',names,'weights',W,'surrogate',...
                        surrogate,'skipchecks',skipchecks,'stream',Stream);
                    end     %parallel loop over i
                end
            else
                parfor i = 1:nTrees
                    sampleidx = 1:length(Y);
                    ibidx = randsample(sampleidx,nboot,SampleWithReplacement);
                    oobidx{i} = setdiff(sampleidx,ibidx);
                    Tree{i} = classregtree2(X(ibidx,:),Y(ibidx,:),...
                        'priorprob',Prior,'cost',Cost,'splitcriterion',...
                        Criterion,'splitmin',splitmin,'minparent',...
                        minparent,'minleaf',minleaf,'nvartosample',...
                        nvartosample,'mergeleaves',Merge,'categorical',...
                        categ,'prune',Prune,'method',Method,'qetoler',...
                        qetoler,'names',names,'weights',W,'surrogate',...
                        surrogate,'skipchecks',skipchecks,'stream',Stream);
                end     %parallel loop over i
            end

            forest.Trees = Tree;
            forest.oobidx = oobidx;
            forest.nTrees = length(forest.Trees);
            if ~RandomForest
                forest.rpm = promat;
            end
        end     %function rpclassificationforest
        
        function [err,Yhats,varargout] = oobpredict(forest,X,Y,treenum)
            if nargin == 3
                treenum = 'last';
            end
            nrows = size(X,1);
            predmat = NaN(nrows,forest.nTrees);
            predcell = cell(nrows,forest.nTrees);
            OOBIndices = forest.oobidx;
            trees = forest.Trees;
            Labels = forest.classname;
            rpm = forest.rpm;
            if isempty(rpm)
                parfor i = 1:forest.nTrees
                    pred_i = num2cell(NaN(nrows,1));
                    pred_i(OOBIndices{i}) = eval(trees{i},X(OOBIndices{i},:));
                    predcell(:,i) = pred_i;
                end
            else
                parfor i = 1:forest.nTrees
                    pred_i = num2cell(NaN(nrows,1));
                    pred_i(OOBIndices{i}) = eval(trees{i},X(OOBIndices{i},:)*rpm{i});
                    predcell(:,i) = pred_i;
                end
            end
            for j = 1:length(Labels)
                predmat(strcmp(predcell,Labels{j})) = j;
            end
            if strcmp(treenum,'every')
                err = NaN(forest.nTrees,1);
                for i = 1:forest.nTrees
                    ensemblepredictions = mode(predmat(:,1:i),2);
                    missing = isnan(ensemblepredictions);
                    Yhats = Labels(ensemblepredictions(~missing));
                    wrong = ~strcmp(Yhats,Y(~missing));
                    err(i) = mean(wrong);
                end
            else
                ensemblepredictions = mode(predmat,2);
                missing = isnan(ensemblepredictions);
                Yhats = Labels(ensemblepredictions(~missing));
                wrong = ~strcmp(Yhats,Y(~missing));
                err = mean(wrong);         
            end
            if length(unique(Y)) == 2
                pos = num2str(max(str2num(char(Y))));
                neg = num2str(min(str2num(char(Y))));
                %varargout{1} = sum(strcmp(predictions(strcmp(pos,Y)),Y(~missing & strcmp(pos,Y))))/sum(strcmp(pos,Y));  %sensitivity
                %varargout{2} = sum(strcmp(predictions(strcmp(pos,Y)),Y(~missing & strcmp(pos,Y))))/sum(strcmp(pos,predictions));    %ppv
                %varargout{3} = sum(strcmp(predictions(strcmp(neg,Y)),Y(~missing & strcmp(neg,Y))))/sum(strcmp(neg,Y));  %specificity
                %varargout{4} = sum(strcmp(predictions(strcmp(neg,Y)),Y(~missing & strcmp(neg,Y))))/sum(strcmp(neg,predictions));    %npv
                varargout{1} = sum(strcmp(Yhats(strcmp(pos,Y)),Y(~missing & strcmp(pos,Y)))); %tp
                varargout{2} = sum(~strcmp(Yhats(strcmp(pos,Y)),Y(~missing & strcmp(pos,Y))));    %fn
                varargout{3} = sum(strcmp(Yhats(strcmp(neg,Y)),Y(~missing & strcmp(neg,Y)))); %tn
                varargout{4} = sum(~strcmp(Yhats(strcmp(neg,Y)),Y(~missing & strcmp(neg,Y))));    %fp
                
            end
        end     %function oobpredict
        
        function Yhats = predict2(forest,X)
            nrows = size(X,1);
            predmat = NaN(nrows,forest.nTrees);
            predcell = cell(nrows,forest.nTrees);
            trees = forest.Trees;
            Labels = forest.classname;
            rpm = forest.rpm;
            
            if isempty(rpm)
                parfor i = 1:forest.nTrees
                    pred_i = eval(trees{i},X);
                    predcell(:,i) = pred_i;
                end
            else
                parfor i = 1:forest.nTrees
                    pred_i = eval(trees{i},X*rpm{i});
                    predcell(:,i) = pred_i;
                end
            end
            for j = 1:length(Labels)
                predmat(strcmp(predcell,Labels{j})) = j;
            end
            
            ensemblepredictions = mode(predmat,2);
            Yhats = Labels(ensemblepredictions);
        end     %function predict2
        
        function Y = predict(forest,X)
            predmat = NaN(size(X,1),forest.nTrees);
            YTree = cell(size(X,1),forest.nTrees);
            Tree = forest.Trees;
            rpm = forest.rpm;
            if isempty(rpm)
                parfor i = 1:forest.nTrees
                    YTree(:,i) = eval(Tree{i},X);
                end
            else
                parfor i = 1:forest.nTrees
                    YTree(:,i) = eval(Tree{i},X*rpm{i});
                end
            end
            classname = forest.classname;
            for j = 1:length(classname)
                predmat(strcmp(YTree,classname{j})) = j;
            end
            ensemblepredictions = mode(predmat,2);
            missing = isnan(ensemblepredictions);
            Y = classname(ensemblepredictions(~missing));
        end     %function predict
    end     %methods
end     %classdef

%----------------------------------------------------
function M = srpmat(d,k,method,varargin)
    if strcmp(method,'old')
        s = varargin{1};
        M = sparse(vec2mat(randsample([-sqrt(s) 0 sqrt(s)],d*k,true,[1/(2*s) 1-1/s 1/(2*s)]),k));
    else
        M = sparse(d,k);
        R = poissrnd(1,1,k);
        %R(R==0) = 1;
        R(R>d) = d;
        rowidx = arrayfun(@(n) randsample(d,n,false),R,'UniformOutput',false);
        for j = 1:k
            M(rowidx{j},j) = 1;
        end
        %M(:,all(M==0,1)) = [];
        if strcmp(method,'guassian')
            M(M~=0) = randn(nnz(M),1);
        elseif strcmp(method,'uniform')
            syms a x
            eqn = symsum(x^2,-a,a)/(2*a) == d;
            a = solve(eqn,a);
            a = double(a);
            a = round(a(a>0));
            binsample = randsample([-1 1],sum(R),true,[0.5 0.5]);
            M(M==1) = binsample;
            M(M==1) = unidrnd(a,sum(M(:)==1),1);
            M(M==-1) = -unidrnd(a,sum(M(:)==-1),1);
        else
            binsample = randsample([-1 1],sum(R),true,[0.5 0.5]);
            M(M==1) = binsample;
        end
    end
end