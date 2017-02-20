close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1)

ps = [2,4,6];
ns{1} = [20,200,400];
ns{2} = [80,400,4000];
ns{3} = [400,2000,4000];
ntest = 10000;
ntrials = 10;

for j = 1:length(ps)
    p = ps(j);
    nOrthants = 2^p;
    nNodes = nOrthants*2 - 1;
    nLevels = p + 1;
    nClasses = nOrthants;
    Tree.Node = (1:nNodes)';
    Tree.CutPoint = NaN(nNodes,1);
    Tree.CutVar = NaN(nNodes,1);
    Tree.Children = NaN(nNodes,2);
    Tree.Parent = NaN(nNodes,1);
    for Level = 1:nLevels
        Nodes = 2^(Level-1):2^Level - 1;
        if Level ~= nLevels
            for k = 1:length(Nodes)
                Node = Nodes(k);
                Tree.CutPoint(Node) = 0;
                Tree.CutVar(Node) = Level;
                Tree.Children(Node,:) = [Node*2,Node*2+1];
                Tree.Parent(Tree.Children(Node,:)) = Node;
            end
        end
    end
    for i = 1:length(ns{j})
        % sample training data
        ntrain = ns{j}(i);
        for Trial = 1:ntrials
            X = rand(ntrain,p)*2 - 1;
            Y = zeros(ntrain,1);
            NodeObservations = cell(nNodes,1);
            NodeObservations{1} = 1:ntrain;
            for Level = 1:nLevels
                Nodes = 2^(Level-1):2^Level - 1;
                if Level ~= nLevels
                    for k = 1:length(Nodes)
                        Node = Nodes(k);
                        Xnode = X(NodeObservations{Node},:);
                        MoveLeft = Xnode(:,Tree.CutVar(Node)) <= Tree.CutPoint(Node);
                        NodeObservations{Tree.Children(Node,1)} = NodeObservations{Node}(MoveLeft);
                        NodeObservations{Tree.Children(Node,2)} = NodeObservations{Node}(~MoveLeft);
                    end
                else
                    for k = 1:length(Nodes)
                        Node = Nodes(k);
                        Y(NodeObservations{Node}) = k;
                    end
                end
            end
           dlmwrite(sprintf('%sRandomerForest/Data/Orthant/Orthant_train_p%d_n%d_trial%d.dat',rerfPath,p,ntrain,Trial),[X,Y])
        end
    end
    % sample test data
    X = rand(ntest,p)*2 - 1;
    Y = zeros(ntest,1);
    Posteriors = zeros(ntest,nClasses);
    NodeObservations = cell(nNodes,1);
    NodeObservations{1} = 1:ntest;
    for Level = 1:nLevels
        Nodes = 2^(Level-1):2^Level - 1;
        if Level ~= nLevels
            for k = 1:length(Nodes)
                Node = Nodes(k);
                Xnode = X(NodeObservations{Node},:);
                MoveLeft = Xnode(:,Tree.CutVar(Node)) <= Tree.CutPoint(Node);
                NodeObservations{Tree.Children(Node,1)} = NodeObservations{Node}(MoveLeft);
                NodeObservations{Tree.Children(Node,2)} = NodeObservations{Node}(~MoveLeft);
            end
        else
            for k = 1:length(Nodes)
                Node = Nodes(k);
                Y(NodeObservations{Node}) = k;
            end
        end
    end
    for c = 1:nClasses
        Posteriors(:,c) = double(Y==c);
    end
   dlmwrite(sprintf('%sRandomerForest/Data/Orthant/Orthant_test_p%d.dat',rerfPath,p),[X,Y])
   dlmwrite(sprintf('%sRandomerForest/Data/Orthant/Orthant_test_posteriors_p%d.dat',rerfPath,p),Posteriors)
end