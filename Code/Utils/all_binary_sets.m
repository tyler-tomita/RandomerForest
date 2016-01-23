function s=all_binary_sets(n)
%
% INPUT: n integer >=0
% OUTPUT: binary coding of all subsets of {1,...,n}

if n==0
    s=zeros(1,0,'int8');
else
    snminus1=all_binary_sets(n-1);
    c=class(snminus1);
    j=size(snminus1,1);
    s=[zeros(j,1,c) snminus1; ...
       ones(j,1,c) snminus1]; 
end