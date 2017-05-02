classdef rpclassificationforest
    
    properties
        oobidx = {}; %indices of out of bag samples for each tree
        Tree = {};  %classification trees contained in the ensemble
        nTrees = []; %number of trees in the ensemble
        classname;
        ForestMethod;
        RandomMatrix;
        Rescale;
%         NumVars = [];
        priors = [];
        rotmat = [];
        RotVars = logical([]);
        rpm = [];
        rho = [];
    end
    
    methods
        function forest = rpclassificationforest(X,Y,varargin)
            %class contstructor for RandomerForest object
            %
            %nTrees: number of trees in ensemble
            %
            %X: n x d matrix where n is number of samples and d is
            %
            %number of dimensions (predictor variables)
            %
            %Y: n x 1 cell string of class labels
            %
            %Optional Arguments:
                %nvartosample: if 'RandomForest' is true, then this is the
                %number of variables subsampled when splitting. Otherwise,
                %this is the dimension of the subspace randomly projected
                %into
                %
                %s: s a parameter that specifies the sparsity of the random
                %projection matrix. Sparsity is computed as 1/(2*s). Only
                %used if RandomMatrix is set to 'dense'
                %
                %mdiff: string 'all', 'node' or 'off'. Specifying 'all' or
                %'node' allows the full-sample or node-sample
                %class-conditional differences in means to be sampled as
                %projections
                %
                %RandomMatrix: string specifying the method for sampling
                %the random projection matrix. Options are 'dense' (dense
                %nonzeros sampled from [-1,1]), 'sparse' (nonzeros
                %sampled from {-1,1}), and 'frc' (Breiman's Forest-RC). For
                %dense, sparisty is controlled by the parameter 's'.
                %
                %RandomForest: logical true or false (default). Setting to
                %true performs traditional random forest
                %
                %Rescale: logical true or false (default). Setting to true
                %passes the data to marginal ranks prior to any computing
                %
                %NWorkers: number of workers for parallel computing
                %(default = 1)
                %
                %rotate: logical true or false (default). Setting to true
                %uniformly randomly rotates the data for each tree prior to
                %fitting
                %
                %p: probability of sampling each of the K-1 mean difference
                %vectors, where K is the number of classes
                %
                %dx: size of expanded feature set in one random projection
                %per tree version
            %
            %Example:
            %
            %Train a RerF classifier consisting of 500 trees using dense random projections with
            %sparsity = 2/3 (s = 3) and allowing mean difference
            %projections to be sampled. The 'node' option computes sample
            %means using only samples at the current split node. Connect to
            %6 workers, which allows 6 trees to be constructed in parallel.
            %Compute out of bag error achieved at the final ('last') tree.
            %
            %rerf = rpclassificationforest(500,X,Y,'s',3,'mdiff','node','NWorkers',6);
            %
            %err = oobpredict(rerf,X,Y,'last');
                
            if ~iscell(Y)
                Y = cellstr(num2str(Y));
            end
            forest.classname = unique(Y);
            forest = growTrees(forest,X,Y,varargin{:});
        end     %class constructor
        
        function forest = growTrees(forest,X,Y,varargin)
            okargs =   {'priorprob' 'cost'    'splitcriterion'  'splitmin'...
                        'minparent' 'minleaf'   'nvartosample'...
                        'mergeleaves'   'categorical' 'prune' 'method' ...
                        'qetoler'   'names'   'weights' 'surrogate'...
                        'skipchecks'    'stream'    'fboot'...
                        'SampleWithReplacement' 'rho' 'mdiff' 'ForestMethod'...
                        'RandomMatrix'  'Rescale'   'NWorkers'  'Stratified'...
                        'nmix'  'rotate'    'p' 'dprime'    'nTrees'    'dx'...
                        'AdjustmentFactors'     'DownsampleNode'    'MaxNodeSize'};
            defaults = {[]  []  'gdi'   []  2  1   ceil(size(X,2)^(2/3))...
                        'off'    []  'off'    'classification'  1e-6    {}...
                        []  'off'   false  []  1    true   1/size(X,2)    'off'...
                        'rerf'  'sparse'    'off' 1   true   2   false  []...
                        []  500   size(X,2) []  false   100};
            [Prior,Cost,Criterion,splitmin,minparent,minleaf,...
                nvartosample,Merge,categ,Prune,Method,qetoler,names,W,...
                surrogate,skipchecks,Stream,fboot,...
                SampleWithReplacement,rho,mdiff,ForestMethod,RandomMatrix,...
                Rescale,NWorkers,Stratified,nmix,rotate,p,dprime,nTrees,dx,...
                AdjustmentFactors,DownsampleNode,MaxNodeSize,~,extra] = internal.stats.parseArgs(okargs,defaults,varargin{:});
            
            % If ForestMethod is F-RC and nmix = 1, then just do RF instead
            % since it's the same and faster
            if strcmp(ForestMethod,'rerf') && strcmp(RandomMatrix,'frc')
                if nmix==1
                    ForestMethod = 'rf';
                end
            end
            
            %Convert to double if not already
            if ~isa(X,'double')
                X = double(X);
            end
            
            if ~strcmp(Rescale,'off')
                X = rescale(X,[],Rescale);
                forest.Rescale = Rescale;
            else
                forest.Rescale = 'off';
            end
            
            [n,d] = size(X);
            
            %Check sparsity
            if ischar(rho)
                if ~strcmp(rho,'random')
                    error('Parameter rho must either be ''random'' or type numeric');
                else
                    rho = randi(3,1,nTrees)/d;
                end
            else
                if rho < 1/d
                    rho = repmat(1/d,1,nTrees);
                elseif rho > 1
                    rho = ones(1,nTrees);
                else
                    rho = repmat(rho,1,nTrees);
                end
            end
            
            nclasses = length(forest.classname);
            priors = NaN(1,nclasses);
            for c = 1:nclasses
                priors(c) = sum(strcmp(Y,forest.classname(c)))/length(Y);
            end
            nboot = ceil(fboot*length(Y));
            Tree = cell(nTrees,1);
            oobidx = cell(nTrees,1);
            sampleidx = 1:length(Y);
            
            if d<=500
                RR = zeros(d,d,nTrees);
            else
                RR = zeros(500,500,nTrees);
            end
            RotVars = false(nTrees,d);
            
            % one random projection per tree implementation
%             load Random_matrix_adjustment_factor
            slope = [0.7205 0.7890 0.8143 0.8298 0.8442 0.8600 0.8794 0.8916 0.8922];
            dims = [2 5 10 25 50 100 250 500 1000];
            RM = cell(nTrees,1);
%                 if d <= 10
%                     dx = 2^d;
%                 elseif d > 10 && d <= 100
%                     dx = d^2;
%                 elseif d > 100 && d <= 1000
%                     dx = ceil(d^1.5);
%                 end
            
%             poolobj = gcp('nocreate');
%             if isempty(poolobj);
%                 parpool('local',NWorkers,'IdleTimeout',360);
%             end
            
            %Reduce computational load for matrices with many zero elements
            if ~strcmp(ForestMethod,'rf') && ~rotate && ~issparse(X) && nnz(X)/numel(X) <= 0.01
                X = sparse(X);
            end
            
%             for i = 1:nTrees
            parfor i = 1:nTrees
                %Rotate data?
                if rotate
                    %is d > 500? if so only rotate a subset of 500 of the
                    %dimensions
                    if d<=500
                        RR(:,:,i) = random_rotation(d);
                        Xtree = X*RR(:,:,i);
                    else
                        RR(:,:,i) = random_rotation(500);
                        RotVars(i,:) = ismember(1:d,randperm(d,500));
                        Xtree = X;
                        Xtree(:,RotVars(i,:)) = Xtree(:,RotVars(i,:))*RR(:,:,i);
                    end
                else
                    Xtree = X;
                end
                
                go = true;
                if Stratified
                    while go
                        ibidx = [];
                        for c = 1:nclasses
                            idx = find(strcmp(forest.classname{c},Y));
                            if length(idx) > 1
                                ibidx = cat(2,ibidx,transpose(randsample(idx,ceil(fboot*length(idx)),SampleWithReplacement)));
                            else
                                ibidx(end+1) = idx;
                            end
                        end
                        oobidx{i} = setdiff(sampleidx,ibidx);
                        go = isempty(oobidx{i});
                    end
                else
                    while go
                        ibidx = randsample(sampleidx,nboot,SampleWithReplacement);
                        oobidx{i} = setdiff(sampleidx,ibidx);
                        go = isempty(oobidx{i});
                    end
                end
                
                if ~strcmp(ForestMethod,'rf') && ~strcmp(ForestMethod,'rerf2')
                    Tree{i} = rpclassregtree(Xtree(ibidx,:),Y(ibidx,:),...
                        'priorprob',Prior,'cost',Cost,'splitcriterion',...
                        Criterion,'splitmin',splitmin,'minparent',...
                        minparent,'minleaf',minleaf,'nvartosample',...
                        nvartosample,'mergeleaves',Merge,'categorical',...
                        categ,'prune',Prune,'method',Method,'qetoler',...
                        qetoler,'names',names,'weights',W,'surrogate',...
                        surrogate,'skipchecks',skipchecks,'stream',Stream,...
                        'rho',rho(i),'mdiff',mdiff,'RandomMatrix',RandomMatrix,...
                        'nmix',nmix,'p',p,'dprime',dprime,'DownsampleNode',...
                        DownsampleNode,'MaxNodeSize',MaxNodeSize);
                elseif strcmp(ForestMethod,'rerf2')
                    RM{i} = randmat(d,dx,RandomMatrix,rho(i),nmix,...
                        ceil(dx^(1/interp1(AdjustmentFactors.dims,AdjustmentFactors.slope,d))));
                    Tree{i} = classregtree2(Xtree(ibidx,:)*RM{i},Y(ibidx,:),...
                        'priorprob',Prior,'cost',Cost,'splitcriterion',...
                        Criterion,'splitmin',splitmin,'minparent',...
                        minparent,'minleaf',minleaf,'nvartosample',...
                        min(nvartosample,size(RM{i},2)),'mergeleaves',...
                        Merge,'categorical',categ,'prune',Prune,'method',...
                        Method,'qetoler',qetoler,'names',names,'weights',...
                        W,'surrogate',surrogate,'skipchecks',skipchecks,...
                        'stream',Stream,'DownsampleNode',DownsampleNode,...
                        'MaxNodeSize',MaxNodeSize);
                else
                    Tree{i} = classregtree2(Xtree(ibidx,:),Y(ibidx,:),...
                        'priorprob',Prior,'cost',Cost,'splitcriterion',...
                        Criterion,'splitmin',splitmin,'minparent',...
                        minparent,'minleaf',minleaf,'nvartosample',...
                        nvartosample,'mergeleaves',Merge,'categorical',...
                        categ,'prune',Prune,'method',Method,'qetoler',...
                        qetoler,'names',names,'weights',W,'surrogate',...
                        surrogate,'skipchecks',skipchecks,'stream',Stream,...
                        'DownsampleNode',DownsampleNode,'MaxNodeSize',MaxNodeSize);
                end  
            end     %parallel loop over i
            
            %Compute interpretability as total number of variables split on
%             NumVars = NaN(1,nTrees);
%             if RandomForest
%                 for i = 1:nTrees
%                     NumVars(i) = sum(Tree{i}.var~=0);
%                 end
%             else
%                 for i = 1:nTrees
%                     internalnodes = transpose(Tree{i}.node(Tree{i}.var ~= 0));
%                     TreeVars = zeros(1,length(Tree{i}.node));
%                     for nd = internalnodes
%                         if ~Tree{i}.isdelta(nd)
%                             TreeVars(nd) = nnz(Tree{i}.rpm{nd});
%                         end
%                     end
%                     NumVars(i) = sum(TreeVars);
%                 end
%             end                        
            forest.Tree = Tree;
            forest.oobidx = oobidx;
            forest.nTrees = length(forest.Tree);
            forest.ForestMethod = ForestMethod;
            forest.RandomMatrix = RandomMatrix;
%             forest.NumVars = NumVars;
            forest.priors = priors;
            forest.rho = rho;
            if rotate
                forest.rotmat = RR;
                if d>500
                    forest.RotVars = RotVars;
                end
            end
            if strcmp(ForestMethod,'rerf2')
                forest.rpm = RM;
            end
        end     %function rpclassificationforest
        
%         oobpredict is deprecated
%         function [predcell,err] = oobpredict(forest,X,Y)
%             
%             %Convert to double if not already
%             if ~isa(X,'double')
%                 X = double(X);
%             end
%             
%             if ~strcmp(forest.Rescale,'off')
%                 X = rescale(X,[],forest.Rescale);
%             end
%             [nrows,d] = size(X);
%             predcell = cell(nrows,forest.nTrees);
%             err = NaN(1,forest.nTrees);
%             OOBIndices = forest.oobidx;
%             trees = forest.Tree;
%             rotate = ~isempty(forest.rotmat);
%             if ~strcmp(forest.ForestMethod,'rf')
%                 parfor i = 1:forest.nTrees
%                     if rotate
%                         if d<=500
%                             Xtree = X*forest.rotmat(:,:,i);
%                         else
%                             Xtree = X;
%                             Xtree(:,forest.RotVars(i,:)) = X(:,forest.RotVars(i,:))*forest.rotmat(:,:,i);       
%                         end
%                     else
%                         Xtree = X;
%                     end
%                     pred_i = cell(nrows,1);
%                     pred_i(OOBIndices{i}) = rptreepredict(trees{i},Xtree(OOBIndices{i},:));
%                     predcell(:,i) = pred_i;
%                     err(i) = sum(~strcmp(pred_i(OOBIndices{i}),Y(OOBIndices{i})))/length(OOBIndices{i});
%                 end
%             else
%                 parfor i = 1:forest.nTrees
%                     if rotate
%                         if d<=500
%                             Xtree = X*forest.rotmat(:,:,i);
%                         else
%                             Xtree = X;
%                             Xtree(:,forest.RotVars(i,:)) = Xtree(:,forest.RotVars(i,:))*forest.rotmat(:,:,i);       
%                         end
%                     else
%                         Xtree = X;
%                     end
%                     pred_i = cell(nrows,1);
%                     pred_i(OOBIndices{i}) = eval(trees{i},Xtree(OOBIndices{i},:));
%                     predcell(:,i) = pred_i;
%                     err(i) = sum(~strcmp(pred_i(OOBIndices{i}),Y(OOBIndices{i})))/length(OOBIndices{i});
%                 end
%             end
%         end     %function oobpredict
        
        function scores = rerf_oob_classprob(forest,Xtrain,treenum)
            if nargin == 2
                treenum = 'last';
            end
            
            %Convert to double if not already
            if ~isa(Xtrain,'double')
                Xtrain = double(Xtrain);
            end
            
            if ~strcmp(forest.Rescale,'off')
                Xtrain = rescale(Xtrain,[],forest.Rescale);
            end
            [nrows,d] = size(Xtrain);
            
            Labels = forest.classname;
            nclasses = length(Labels);
            scoremat = NaN(nrows,nclasses,forest.nTrees);
            OOBIndices = forest.oobidx;
            trees = forest.Tree;
            rotate = ~isempty(forest.rotmat);
            RR = forest.rotmat;
            RotVars = forest.RotVars;
            RM = forest.rpm;
            priors = forest.priors;
            if strcmp(forest.ForestMethod,'rerf')
                parfor i = 1:forest.nTrees
                    Xtree = Xtrain(OOBIndices{i},:);
                    score_i = NaN(nrows,nclasses);
                    score_i(OOBIndices{i},:) = rpclassprob(trees{i},Xtree)
                    scoremat(:,:,i) = score_i;
                end
            elseif strcmp(forest.ForestMethod,'rerf2')
                parfor i = 1:forest.nTrees
                    Xtree = Xtrain(OOBIndices{i},:)*RM{i};
                    score_i = NaN(nrows,nclasses);
                    score_i(OOBIndices{i},:) = rfclassprob(trees{i},Xtree)
                    scoremat(:,:,i) = score_i;
                end
            else
                if rotate
                    if d<=500
                        parfor i = 1:forest.nTrees
                            Xtree = Xtrain(OOBIndices{i},:)*RR(:,:,i);
                            score_i = NaN(nrows,nclasses);
                            score_i(OOBIndices{i},:) = rfclassprob(trees{i},Xtree);
                            scoremat(:,:,i) = score_i;
                        end
                    else
                        parfor i = 1:forest.nTrees
                            Xtree = Xtrain(OOBIndices{i},:);
                            Xtree(:,RotVars(i,:)) = Xtree(:,RotVars(i,:))*RR(:,:,i);   
                            score_i = NaN(nrows,nclasses);
                            score_i(OOBIndices{i},:) = rfclassprob(trees{i},Xtree);
                            scoremat(:,:,i) = score_i;
                        end
                    end
                else
                    parfor i = 1:forest.nTrees
                        Xtree = Xtrain(OOBIndices{i},:);
                        score_i = NaN(nrows,nclasses);
                        score_i(OOBIndices{i},:) = rfclassprob(trees{i},Xtree);
                        scoremat(:,:,i) = score_i;
                    end
                end
            end
            if strcmp(treenum,'every')
                scores = NaN(size(scoremat));
                parfor i = 1:forest.nTrees
                    score_i = nanmean(scoremat(:,:,1:i),3);
                    missing = any(isnan(score_i),2);
                    %fprintf('%d\n',size(score_i(missing,:)))
                    %fprintf('%d\n',size(repmat(forest.priors,length(missing),1)))
                    score_i(missing,:) = repmat(priors,sum(missing),1);
                    scores(:,:,i) = score_i;
                end
            else
                scores = nanmean(scoremat,3);
                missing = any(isnan(scores),2);
                scores(missing,:) = repmat(priors,sum(missing),1);
            end
        end     %function rerf_oob_classprob
        
        function scores = rerf_classprob(forest,Xtest,treenum,varargin)
            if nargin == 2
                treenum = 'last';
            end
            
            if nargin == 4;
                Xtrain = varargin{1};
                if ~isa(Xtrain,'double')
                    Xtrain = double(Xtrain);
                end
            end
            
            if ~strcmp(forest.Rescale,'off')
                if nargin < 4
                    error('Training data is required as third input argument for predicting')
                end
                Xtest = rescale(Xtrain,Xtest,forest.Rescale);
            end
            
            %Convert to double if not already
            if ~isa(Xtest,'double')
                Xtest = double(Xtest);
            end
            
            [nrows,d] = size(Xtest);
            
            Labels = forest.classname;
            nclasses = length(Labels);
            scoremat = NaN(nrows,nclasses,forest.nTrees);
            trees = forest.Tree;
            rotate = ~isempty(forest.rotmat);
            RR = forest.rotmat;
            RotVars = forest.RotVars;
            RM = forest.rpm;
            if ~strcmp(forest.ForestMethod,'rf') && ~strcmp(forest.ForestMethod,'rerf2')
                parfor i = 1:forest.nTrees
                    score_i = rpclassprob(trees{i},Xtest)
                    scoremat(:,:,i) = score_i;
                end
            elseif strcmp(forest.ForestMethod,'rerf2')
                parfor i = 1:forest.nTrees
                    Xtree = Xtest*RM{i};
                    score_i = rfclassprob(trees{i},Xtree);
                    scoremat(:,:,i) = score_i;
                end
            else
                if rotate
                    if d<=500
                        parfor i = 1:forest.nTrees
                            Xtree = Xtest*RR(:,:,i);
                            score_i = rfclassprob(trees{i},Xtree);
                            scoremat(:,:,i) = score_i;
                        end
                    else
                        parfor i = 1:forest.nTrees
                            Xtree = Xtest;
                            Xtree(:,RotVars(i,:)) = Xtree(:,RotVars(i,:))*RR(:,:,i);       
                            score_i = rfclassprob(trees{i},Xtree);
                            scoremat(:,:,i) = score_i;
                        end
                    end
                else
                    parfor i = 1:forest.nTrees
                        Xtree = Xtest;
                        score_i = rfclassprob(trees{i},Xtree);
                        scoremat(:,:,i) = score_i;
                    end
                end
            end
            if strcmp(treenum,'every')
                scores = NaN(size(scoremat));
                parfor i = 1:forest.nTrees
                    score_i = mean(scoremat(:,:,1:i),3);
                    scores(:,:,i) = score_i;
                end
            elseif strcmp(treenum,'individual')
                scores = scoremat;
            else
                scores = mean(scoremat,3);
            end
        end     %function rerf_classprob
        
%         DEPRECATED 
%         function Y = predict(forest,Xtest,varargin)
%                         
%             %Convert to double if not already
%             if ~isa(Xtest,'double')
%                 Xtest = double(Xtest);
%             end
%             
%             if nargin == 3;
%                 Xtrain = varargin{1};
%                 if ~isa(Xtrain,'double')
%                     Xtrain = double(Xtrain);
%                 end
%             end
%             
%             if ~strcmp(forest.Rescale,'off')
%                 if nargin < 3
%                     error('Training data is required as third input argument for predicting')
%                 end
%                 Xtest = rescale(Xtrain,Xtest,forest.Rescale);
%             end
%             
%             [n,d] = size(Xtest);
%             predmat = NaN(n,forest.nTrees);
%             YTree = cell(n,forest.nTrees);
%             Tree = forest.Tree;
%             rotate = ~isempty(forest.rotmat);
%             
%             if ~strcmp(forest.ForestMethod,'rf')
%                 parfor i = 1:forest.nTrees
%                     if rotate
%                         if d<=500
%                             Xtree = Xtest*forest.rotmat(:,:,i);
%                         else
%                             Xtree = Xtest;
%                             Xtree(:,forest.RotVars(i,:)) = Xtree(:,forest.RotVars(i,:))*forest.rotmat(:,:,i);       
%                         end
%                     else
%                         Xtree = Xtest;
%                     end
%                     YTree(:,i) = rptreepredict(Tree{i},Xtree);
%                 end
%             else
%                 parfor i = 1:forest.nTrees
%                     if rotate
%                         if d<=500
%                             Xtree = Xtest*forest.rotmat(:,:,i);
%                         else
%                             Xtree = Xtest;
%                             Xtree(:,forest.RotVars(i,:)) = Xtree(:,forest.RotVars(i,:))*forest.rotmat(:,:,i);       
%                         end
%                     else
%                         Xtree = Xtest;
%                     end
%                     YTree(:,i) = eval(Tree{i},Xtree);
%                 end
%             end
%             Labels = forest.classname;
%             for j = 1:length(Labels)
%                 predmat(strcmp(YTree,Labels{j})) = j;
%             end
%             ensemblepredictions = mode(predmat,2);
%             missing = isnan(ensemblepredictions);
%             Y = Labels(ensemblepredictions(~missing));
%         end     %function predict
        
%             DEPRECATED
%             function sp = db_sparsity(forest)
%             %sparsity of decision boundary computed as sum #variables used
%             %over all nodes
%             
%             sp = 0;
%             for i = 1:forest.nTrees
%                 Tree = forest.Tree{i};
%                 if ~strcmp(forest.ForestMethod,'rf')
%                     internalnodes = Tree.node(Tree.var~=0);
%                     for node = internalnodes'
%                         sp = sp + sum(Tree.rpm{node}~=0);
%                     end
%                 else
%                     sp = sp + sum(Tree.var~=0);
%                 end
%             end
%         end
        
        function Depth = forest_depth(forest)
            Depth = NaN(forest.nTrees,1);
            trees = forest.Tree;
            parfor i = 1:forest.nTrees
                Depth(i) = tree_depth(trees{i});
            end
        end % forest_depth
        
        function tree_strength(forest,Xtest,Ytest,varargin)
            
            if nargin == 4;
                Xtrain = varargin{1};
                if ~isa(Xtrain,'double')
                    Xtrain = double(Xtrain);
                end
            end
            
            if ~strcmp(forest.Rescale,'off')
                if nargin < 4
                    error('Training data is required as third input argument for predicting')
                end
                Xtest = rescale(Xtrain,Xtest,forest.Rescale);
            end
            
            %Convert to double if not already
            if ~isa(Xtest,'double')
                Xtest = double(Xtest);
            end
            
            [nrows,d] = size(Xtest);
            
            Labels = forest.classname;
            nclasses = length(Labels);
            scoremat = NaN(nrows,nclasses,forest.nTrees);
            trees = forest.Tree;
            rotate = ~isempty(forest.rotmat);
            RR = forest.rotmat;
            RotVars = forest.RotVars;
            RM = forest.rpm;
            if ~strcmp(forest.ForestMethod,'rf') && ~strcmp(forest.ForestMethod,'rerf2')
                parfor i = 1:forest.nTrees
                    score_i = rpclassprob(trees{i},Xtest)
                    scoremat(:,:,i) = score_i;
                end
            elseif strcmp(forest.ForestMethod,'rerf2')
                parfor i = 1:forest.nTrees
                    Xtree = Xtest*RM{i};
                    score_i = rfclassprob(trees{i},Xtree);
                    scoremat(:,:,i) = score_i;
                end
            else
                if rotate
                    if d<=500
                        parfor i = 1:forest.nTrees
                            Xtree = Xtest*RR(:,:,i);
                            score_i = rfclassprob(trees{i},Xtree);
                            scoremat(:,:,i) = score_i;
                        end
                    else
                        parfor i = 1:forest.nTrees
                            Xtree = Xtest;
                            Xtree(:,RotVars(i,:)) = Xtree(:,RotVars(i,:))*RR(:,:,i);       
                            score_i = rfclassprob(trees{i},Xtree);
                            scoremat(:,:,i) = score_i;
                        end
                    end
                else
                    parfor i = 1:forest.nTrees
                        Xtree = Xtest;
                        score_i = rfclassprob(trees{i},Xtree);
                        scoremat(:,:,i) = score_i;
                    end
                end
            end
            if strcmp(treenum,'every')
                scores = NaN(size(scoremat));
                parfor i = 1:forest.nTrees
                    score_i = mean(scoremat(:,:,1:i),3);
                    scores(:,:,i) = score_i;
                end
            else
                scores = mean(scoremat,3);
            end
        end     %tree_strength
    end     %methods
end     %classdef
