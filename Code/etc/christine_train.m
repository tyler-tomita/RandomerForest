%AutoML Challenge christine data

close all
clear
clc

fpath = mfilename('fullpath');
findex = strfind(fpath,'/');
rootDir=fpath(1:findex(end-1));
p = genpath(rootDir);
gits=strfind(p,'.git');
colons=strfind(p,':');
for i=0:length(gits)-1
endGit=find(colons>gits(end-i),1);
p(colons(endGit-1):colons(endGit)-1)=[];
end
addpath(p);

X = dlmread('~/Documents/MATLAB/christine/christine_train.data');
Y = cellstr(num2str(dlmread('~/Documents/MATLAB/christine/christine_train.solution')));

d = size(X,2);
ntrees = 1000;
embed = round(d.^[1/4 1/3 1/2 2/3 3/4]);

rferr = NaN(1,length(embed));
f1err = NaN(1,length(embed));
f2err = NaN(1,length(embed));
f3err = NaN(1,length(embed));
f4err = NaN(1,length(embed));
trf = NaN(1,length(embed));
tf1 = NaN(1,length(embed));
tf2 = NaN(1,length(embed));
tf3 = NaN(1,length(embed));
tf4 = NaN(1,length(embed));

parpool;
for i = 1:length(embed)
    nvartosample = embed(i);
    
    tic
    rf = rpclassificationforest(ntrees,X,Y,'nvartosample',nvartosample,'RandomForest',true);
    trf(i) = toc;
    rferr(i) = oobpredict(rf,X,Y,'last');
    clear rf


    fprintf('Random Forest complete\n')

    tic
    f1 = rpclassificationforest(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','old');
    tf1(i) = toc;
    f1err(i) = oobpredict(f1,X,Y,'last');
    clear f1

    fprintf('TylerForest complete\n')

    tic
    f2 = rpclassificationforest(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','new');
    tf2(i) = toc;
    f2err(i) = oobpredict(f2,X,Y,'last');
    clear f2

    fprintf('TylerForest+ complete\n')

    tic
    f3 = rpclassificationforest(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
    tf3(i) = toc;
    f3err(i) = oobpredict(f3,X,Y,'last');
    clear f3

    fprintf('TylerForest+meandiff complete\n')
    
    tic
    f4 = rpclassificationforest(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new','Robust',true);
    tf4(i) = toc;
    f4err(i) = oobpredict(f4,X,Y,'last');
    clear f4

    fprintf('Robust TylerForest+meandiff complete\n')
end

plot(embed,rferr,'-bs',embed,f1err,'-rs',embed,f2err,'-gs',embed,f3err,'-cs',embed,f4err,'-ms')
xlabel('Number of Embedded Dimensions')
ylabel('OOB Error')
legend('Random Forest','Dense Randomer Forest','Sparse Randomer Forest','Sparse Randomer Forest Mean Diff','Robust Sparse Randomer Forest Mean Diff')
fname = 'christine_train_lhat';
save_fig(gcf,fname);
plot(embed,trf,'-bs',embed,tf1,'-rs',embed,tf2,'-gs',embed,tf3,'-cs',embed,tf4,'-ms')
xlabel('Number of Embedded Dimensions')
ylabel('Training Time (sec)')
legend('Random Forest','Dense Randomer Forest','Sparse Randomer Forest','Sparse Randomer Forest Mean Diff','Robust Sparse Randomer Forest Mean Diff')
fname = 'christine_train_time';
save_fig(gcf,fname)