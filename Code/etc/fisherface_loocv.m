function [Lhat,linclass] = fisherface_loocv(X,Y)
    [~,pc] = pca(X);
    linclass = fitcdiscr(pc,Y,'DiscrimType','Linear','FillCoeffs','off','CrossVal','on','Leaveout','on');
    Lhat = kfoldLoss(linclass);
end