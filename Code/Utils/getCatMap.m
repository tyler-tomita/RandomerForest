function iFeatureNum = getCatMap(filename,X)
% GETFOLDS reads a csv file in which each line enumerates the column
% indices of an n-by-p data matrix corresponding to the same categorical
% feature after binary expansion. Used for CCF.

[n,p] = size(X);

iFeatureNum = 1:p;

fid = fopen(filename);

i = 0;
while ~feof(fid)
    i = i + 1;
    catMap = cellfun(@str2num,strsplit(fgetl(fid), ','));
    if i == 1
        featureIdx = catMap(1);
    end
    iFeatureNum(catMap) = featureIdx;
    featureIdx = featureIdx + 1;
end