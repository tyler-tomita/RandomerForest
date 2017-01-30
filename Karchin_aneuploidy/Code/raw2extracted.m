%% Extract features and class labels from Rachel Karchin's raw aneuploidy data and write to tab-delimited files

close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

Folds = dir('~/Karchin_aneuploidy/Data/raw/*train*.csv');

nFolds = length(Folds);

for i = 1:nFolds
    
    % get fold id from file name
    FoldId = regexp(Folds(i).name,'\w*(?=.train.csv)','match');
    FoldId = FoldId{1};
    
    fprintf('Fold %s (%d of %d)\n',FoldId,i,nFolds)
    
    % load training data and extract relevant columns, placing features in
    % the first p columns and class labels as the last column
    TrainSet = {};
    fid = fopen(['~/Karchin_aneuploidy/Data/raw/' FoldId '.train.csv'],'r');
    TextLine = fgetl(fid);
    while ischar(TextLine)
        SplitString = strsplit(TextLine,',');
        TrainSet = [TrainSet;SplitString([8:85,1])];
        TextLine = fgetl(fid);
    end
    fclose(fid);
    TrainSet = cellfun(@str2num,TrainSet);
    
    % load test data and extract relevant columns, placing features in
    % the first p columns and class labels as the last column
    TestSet = {};
    fid = fopen(['~/Karchin_aneuploidy/Data/raw/' FoldId '.test.csv'],'r');
    TextLine = fgetl(fid);
    while ischar(TextLine)
        SplitString = strsplit(TextLine,',');
        TestSet = [TestSet;SplitString([8:85,1])];
        TextLine = fgetl(fid);
    end
    fclose(fid);
    TestSet = cellfun(@str2num,TestSet);
    
    % write extracted data to tab-delimited 
    dlmwrite(['~/Karchin_aneuploidy/Data/extracted/' FoldId '.train.dat'],...
        TrainSet,'delimiter','\t','precision',12);
    dlmwrite(['~/Karchin_aneuploidy/Data/extracted/' FoldId '.test.dat'],...
        TestSet,'delimiter','\t','precision',12);
end