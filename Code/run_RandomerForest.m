function [Lhat,Time,Cl] = run_RandomerForest(FileName)

    X = dlmread(FileName,'\t',1,1);
    Y = cellstr(num2str(X(:,end)));
    X(:,end) = [];
    d = size(X,2);
    embed = ceil(d.^[1/2 2/3 3/4]);
    ntrees = 1000;
    
    parpool;
    for i = 1:length(embed)
        nvartosample = embed(i);

        tic
        Cl.rf(i) = rpclassificationforest(ntrees,X,Y,'nvartosample',nvartosample,'RandomForest',true);
        Time.rf(i) = toc;
        Lhat.rf(i) = oobpredict(Cl.rf(i),X,Y,'last');

        tic
        Cl.f1(i) = rpclassificationforest(ntrees,X,Y,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','dgaussian');
        Time.f1(i) = toc;
        Lhat.f1(i) = oobpredict(Cl.f1(i),X,Y,'last');

        tic
        Cl.f2(i) = rpclassificationforest(ntrees,X,Y,'nvartosample',nvartosample,'mdiff','off','sparsemethod','new');
        Time.f2(i) = toc;
        Lhat.f2(i) = oobpredict(Cl.f2(i),X,Y,'last');

        tic
        Cl.f3(i) = rpclassificationforest(ntrees,X,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        Time.f3(i) = toc;
        Lhat.f3(i) = oobpredict(Cl.f3(i),X,Y,'last');

        tic
        Cl.f4(i) = rpclassificationforest(ntrees,X,Y,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new','Robust',true);
        Time.f4(i) = toc;
        Lhat.f4(i) = oobpredict(Cl.f4(i),X,Y,'last');
    end
end