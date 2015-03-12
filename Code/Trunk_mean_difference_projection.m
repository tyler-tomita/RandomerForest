%Trunk mean difference projection 

clear
close all
clc

ns = [1000];
dims = [1000];
for i = 1:length(dims)
    d = dims(i);
    d_idx = 1:d;
    mu1 = 1./d_idx;
    mu0 = -1*mu1;
    Mu = cat(1,mu0,mu1);
    Sigma = 1*speye(d);
    Class = [0;1];
    b = .2;
    l = [.1 .3 .5 .7];
    h = .7;
    w = .15;
    h_fig = figure(1);
    %set(h_fig,'Position',[0 0 560 140])
    obj = gmdistribution(Mu,Sigma);
    for j = 1:length(ns)
        n = ns(j);
        [X,idx] = random(obj,n);
        Y = Class(idx);
        md = transpose(mean(X(Y==1,:))) - transpose(mean(X(Y==0,:)));
        Xproj = X*md;
        panel = i + (j-1)*length(dims);
        subplot(length(ns),length(dims),panel)
        plot(Xproj(Y==0),zeros(sum(Y==0),1),'bo',Xproj(Y==1),zeros(sum(Y==1,1)),'rx')
        xlabel('mudiff projection')
        %legend('class 0','class 1')
        ax = gca;
        ax.YTick = [];
     %   ax.Position = [l(i) b w h];
        title(sprintf('n = %0.0f, d = %0.0f',n,d))
    end
end

%set(gcf,'PaperOrientation','landscape')
%set(gcf,'NextPlot','add');
axes;
%h = title('Trunk');
set(gca,'Visible','off');
%set(h,'Visible','on');
fname = '~/LOVEFest/Figures/Trunk_mean_difference_projection_n1000';
save_fig(gcf,fname)