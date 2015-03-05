%Make Fig1 for LOVEFest (Randomer Forest)

%Outline:
%Loop over the three existing figures
%Format the figure properties to make it aesthetically pleasing
%Copy to 1x3 subplot

close all
clear
clc

Colors = linspecer(5,'sequential');
Fig_Color = [1 1 1];
LineWidth = 1.5;
Marker = 'none';
Title = {'Trunk' 'Parity' 'Multimodal'};
Units = 'pixels';
FigPosition = [0 140 1150 650];
Axis_Left = [50 350 650 50 350 650];
Axis_Bottom = [350 350 350 50 50 50];
Axis_Width = 250;
Axis_Height = 250;
Legend_Left = 950;
Legend_Bottom = [475 175];
Legend_Width = 150;
Legend_Height = 100;
MarkerSize = 14;
Box = 'off';

BasePath = '~/LOVEFest/Figures/fig/';
Filename = {'trunk_ooberror_vs_d_n100_var1_embed100_ntrees1000_ntrials10.fig'...
 'Parity_ooberror_vs_d_n100_var1_embed100_ntrees1000_ntrials10.fig'...
 'Multimodal_ooberror_vs_d_n100_var1_ntrees1000_ntrials10.fig'};
BayesFigs = {'Trunk_bayes_error_vs_d.fig' 'Parity_bayes_error.fig' 'Multimodal_bayes_error.fig'};


%Copy bayes error line to the axes containing Lhat for the different
%classifiers
for i = 1:length(Filename)
    h{i} = openfig(strcat(BasePath,Filename{i}),'invisible');
    grid off
    h_err = openfig(BayesFigs{i},'invisible');
    ax_old = get(h_err,'CurrentAxes');
    ax_new = get(h{i},'CurrentAxes');
    copyobj(allchild(ax_old),ax_new);
end

h{4} = figure('Visible','On');
set(h{4},'Position',FigPosition,'Color',Fig_Color)

for i = 1:length(h)-1
    ax_old = get(h{i},'CurrentAxes');
    ax_new = subplot(2,3,i);
    copyobj(allchild(ax_old),ax_new);
    h_lines = allchild(ax_new);
    xmax = zeros(1,5);
    ymax = zeros(1,5);
    for j = 1:5
        set(h_lines(j),'Color',Colors(j,:),'linewidth',LineWidth,'Marker',Marker,'MarkerFaceColor',Colors(j,:),'MarkerEdgeColor',Colors(j,:))
        xmax(j) = max(get(h_lines(j),'XData'));
        ymax(j) = max(get(h_lines(j),'YData'));
    end
    xmax = max(xmax);
    ymax = max(ymax);
    XTick = logspace(0,log10(xmax),log10(xmax)+1);
    set(ax_new,'XLim',[0 xmax],'YLim',[0 ymax],'XScale','log','XTick',XTick,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth,'Units',Units,'Position',[Axis_Left(i) Axis_Bottom(i) Axis_Width Axis_Height])
    title(Title{i})
    xlabel('Number of Ambient Dimensions')
    ylabel('L hat')
    hL = legend(ax_new,'Random Forest','Dense Randomer Forest','Sparse Randomer Forest','Sparse Randomer Forest w/ Mean Diff','Bayes Error');
    legend(ax_new,'hide')
    get(ax_new,'Position');
end
set(hL,'Units',Units,'Position',[Legend_Left Legend_Bottom(1) Legend_Width Legend_Height],'Visible','On','Box','Off')

%Scatter Plots
n = 100;
d = 2;
d_idx = 1:d;
mu1 = 1./d_idx;
mu0 = -1*mu1;
mu = cat(1,mu1,mu0);
Sigma = ones(1,d);
obj = gmdistribution(mu,Sigma);
[X,idx] = random(obj,n);
Y = idx;
Y(idx==2) = 0;
ax = subplot(2,3,4);
plot(X(Y==0,1),X(Y==0,2),'.',X(Y==1,1),X(Y==1,2),'.')
xlabel('X1')
ylabel('X2')
ch = allchild(gca);
set(ch(1),'MarkerSize',MarkerSize,'Color','r')
set(ch(2),'MarkerSize',MarkerSize,'Color','b')
set(gca,'XLim',[-4 4],'YLim',[-3 3],'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth,'Units',Units,'Position',[Axis_Left(4) Axis_Bottom(4) Axis_Width Axis_Height])

X = sparse(n,d);
Sigma = ones(1,d);
nones = randi(d,n,1);
Y = mod(nones,2);
Ystr = cellstr(num2str(Y));
Mu = sparse(n,d);
for j = 1:n
    onesidx = randsample(1:d,nones(j),false);
    Mu(j,onesidx) = 1;
end
for j = 1:n
    X(j,:) = mvnrnd(Mu(j,:),Sigma);
end
ax = subplot(2,3,5);
plot(X(Y==0,1),X(Y==0,2),'.',X(Y==1,1),X(Y==1,2),'.')
xlabel('X1')
ylabel('X2')
ch = allchild(gca);
set(ch(1),'MarkerSize',MarkerSize,'Color','r')
set(ch(2),'MarkerSize',MarkerSize,'Color','b')
set(gca,'XLim',[-4 4],'YLim',[-3 3],'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth,'Units',Units,'Position',[Axis_Left(5) Axis_Bottom(5) Axis_Width Axis_Height])

nclasses = 4;
Classes = 1:nclasses;
ncomponents = 2;    %# of mixture components per class
J = nclasses*ncomponents; %total number of mixture components
p = ones(1,J)/J;    %mixture probabilities
Class = repmat(transpose(1:nclasses),1,ncomponents);
Class = Class(:);
Class = Class(randperm(length(Class)));
df = 10*d;   %degrees of freedom for inverse wishart
nvartosample = ceil(d^(2/3));
Mu = zeros(J,d);
Sigma = zeros(d,d,J);
for j = 1:J
Mu(j,:) = mvnrnd(zeros(1,d),eye(d));
Sigma(:,:,j) = iwishrnd(eye(d),df)*(df-d-1);
end
obj = gmdistribution(Mu,Sigma,p);
[X,idx] = random(obj,n);
Y = Class(idx);
ax = subplot(2,3,6);
plot(X(Y==1,1),X(Y==1,2),'.',X(Y==2,1),X(Y==2,2),'.',X(Y==3,1),X(Y==3,2),'.',X(Y==4,1),X(Y==4,2),'.')
xlabel('X1')
ylabel('X2')
ch = allchild(gca);
set(ch(1),'MarkerSize',MarkerSize,'Color','r')
set(ch(2),'MarkerSize',MarkerSize,'Color','b')
set(ch(3),'MarkerSize',MarkerSize,'Color','g')
set(ch(4),'MarkerSize',MarkerSize,'Color','m')
set(gca,'XLim',[-4 4],'YLim',[-3 3],'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth,'Units',Units,'Position',[Axis_Left(6) Axis_Bottom(6) Axis_Width Axis_Height])
hL = legend('Class 1','Class 2','Class 3','Class 4');
set(hL,'Units',Units,'Position',[Legend_Left Legend_Bottom(2) Legend_Width Legend_Height],'Visible','On','Box','Off')
fname = 'Fig1_Lhat';
save_fig(gcf,fname)