%% Print names of each dataset fold to a file

close all
clear
clc

rng(1);

Folds = dir('~/Karchin_aneuploidy/Data/raw/*train*.csv');

nFolds = length(Folds);

fid = fopen('~/Karchin_aneuploidy/Data/Names.txt','w');

for i = 1:nFolds
    
    % get fold id from file name
    FoldId = regexp(Folds(i).name,'\w*(?=.train.csv)','match');
    FoldId = FoldId{1};
    fprintf(fid,'%s\n',FoldId);
end