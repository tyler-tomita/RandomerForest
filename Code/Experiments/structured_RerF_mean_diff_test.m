clear
close all
clc

X = loadMNISTImages('train-images-idx3-ubyte');
X = sparse(X');
Y = loadMNISTLabels('train-labels-idx1-ubyte');
[n,p] = size(X);
ih = sqrt(p);
iw = ih;

Labels = unique(Y);
K = length(Labels);
pairs = zeros(K-1,2);
npairs = K-1;
pairs(:,1) = 1:npairs;
pairs(:,2) = pairs(:,1) + 1;
mu_diff = zeros(p,npairs);

%compute class conditional difference in means 
for i = 1:npairs
    mu_diff(:,i) = transpose(mean(X(Y==Labels(pairs(i,2)),:)) - mean(X(Y==Labels(pairs(i,1)),:)));
end

[A,ClassPairIdx] = structured_rp(ih,iw,[],[],10000,mu_diff);