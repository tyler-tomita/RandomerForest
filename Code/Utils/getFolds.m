function fold = getFolds(filename)
% GETFOLDS reads a csv file in which each line enumerates the indices of a
% data partition for cross-validation

fold = {};

fid = fopen(filename);

i = 0;
while ~feof(fid)
    i = i + 1;
    fold{i} = cellfun(@str2num,strsplit(fgetl(fid), ','));
end