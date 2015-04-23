%Make Fig1 for LOVEFest (Randomer Forest)

close all
clear
clc

stream = RandStream('mt19937ar','Seed',10);
RandStream.setGlobalStream(stream);

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

Colors = linspecer(5,'sequential');
ScatterColors = {'r' 'b' 'g' 'm'};
Fig_Color = [1 1 1];
LineWidth = 1.5;
Marker = 'none';
Title = {'Trunk' 'Parity' 'Multimodal'};
Units = 'pixels';
%FigPosition = [0 140 1150 1150];
FigPosition = [0 140 1150 650];
left = [50 350 650];
bottom = [350 50];
Axis_Left = repmat(left,1,2);
Axis_Bottom = cat(2,repmat(bottom(1),1,3),repmat(bottom(2),1,3));
Axis_Width = 250;
Axis_Height = 250;
Legend_Width = [150 75];
Legend_Height = [100 round(4/5*100)];
Legend_Left = [950 925];
Legend_Bottom = [round(bottom(1)+Axis_Height/2-Legend_Height(1)/2) round(bottom(2)+Axis_Height/2-Legend_Height(2)/2)];
MarkerSize = 14;
Box = 'off';
Box_Legend = 'off';

Level_curve = 1;

BasePath = '~/LOVEFest/Figures/fig/';
Filename = {'Trunk_ooberror_vs_d_n100_var1_ntrees1000_ntrials10_v2.fig'...
 'Parity_ooberror_vs_d_n100_ntrees1000_ntrials10_v5.fig'...
 'Multimodal_ooberror_vs_d_n100_var1_ntrees1000_ntrials10_v2.fig'};
BayesFigs = {'Trunk_bayes_error_vs_d_v2.fig' 'Parity_bayes_error_v4.fig' 'Multimodal_bayes_error.fig'};


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
set(h{4},'Position',FigPosition,'PaperOrientation','landscape','PaperUnits','inches','PaperSize',[8.5*2 11*2],'PaperPositionMode','auto','Color',Fig_Color)

for i = 1:length(h)-1
    ax_old = get(h{i},'CurrentAxes');
    ax_new = subplot(2,3,i);
    copyobj(allchild(ax_old),ax_new);
    h_lines = allchild(ax_new);
    xmax = zeros(1,5);
    ymax_ln = zeros(1,5);
    for j = 1:5
        set(h_lines(j),'Color',Colors(j,:),'linewidth',LineWidth,'Marker',Marker,'MarkerFaceColor',Colors(j,:),'MarkerEdgeColor',Colors(j,:))
        xmax(j) = max(get(h_lines(j),'XData'));
        ymax_ln(j) = max(get(h_lines(j),'YData'));
    end
    xmax = max(xmax);
    ymax_ax = max(ymax_ln);
    ymax_idx = find(ymax_ln==ymax_ax);
    ymax_ax = ymax_ax + h_lines(ymax_idx).UData(find(h_lines(ymax_idx).YData==ymax_ln(ymax_idx)));
    XTick = logspace(0,log10(xmax),log10(xmax)+1);
    set(ax_new,'XLim',[0 10^(log10(xmax)+0.1)],'YLim',[0 ymax_ax],'XScale','log','XTick',XTick,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth,'Units',Units,'Position',[Axis_Left(i) Axis_Bottom(i) Axis_Width Axis_Height])
    title(Title{i})
    xlabel('Number of Ambient Dimensions')
    ylabel('L hat')
    hL = legend(ax_new,'Random Forest','Dense Randomer Forest','Sparse Randomer Forest','Sparse Randomer Forest w/ Mean Diff','Bayes Error');
    legend(ax_new,'hide')
    get(ax_new,'Position');
end
set(hL,'Units',Units,'Position',[Legend_Left(1) Legend_Bottom(1) Legend_Width(1) Legend_Height(1)],'Visible','On','Box',Box_Legend)

%Scatter Plots
n = 100;
d = 2;
Class = [0 1];
d_idx = 1:d;
mu1 = 1./sqrt(d_idx);
mu0 = -1*mu1;
Mu = cat(1,mu0,mu1);
Sigma = 1*speye(d);
obj = gmdistribution(Mu,Sigma);
[X,idx] = random(obj,n);
Y = Class(idx);
ax = subplot(2,3,4);
for j = 1:length(Class)
    plot(X(Y==Class(j),1),X(Y==Class(j),2),'.',...
        'MarkerSize',MarkerSize,...
        'Color',ScatterColors{j})
	hold on
end
for j = 1:size(Mu,1)
    X_level_curve{j} = bvn_level_curve(Mu(j,:),full(Sigma),Level_curve,200);
    plot(X_level_curve{j}(:,1),X_level_curve{j}(:,2),'--',...
        'Color',ScatterColors{j},...
        'LineWidth',LineWidth)
    hold on
end
xlabel('X1')
ylabel('X2')
set(gca,'XLim',[-4 4],'YLim',[-4 4],'YTick',-4:2:4,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth,'Units',Units,'Position',[Axis_Left(4) Axis_Bottom(4) Axis_Width Axis_Height])

X = sparse(n,d);
Sigma = 1/32*speye(d);
Mu = sparse(n,d);

for j = 1:n
    Mu(j,:) = binornd(1,0.5,1,d);
    X(j,:) = mvnrnd(Mu(j,:),Sigma);
end
nones = sum(Mu,2);
Y = mod(nones,2);
ax = subplot(2,3,5);
for j = 1:length(Class)
    plot(X(Y==Class(j),1),X(Y==Class(j),2),'.',...
        'MarkerSize',MarkerSize,...
        'Color',ScatterColors{j})
	hold on
end
Mu = {[0 0] [1 1] [1 0] [0 1]};
Class = [1 1 2 2];
for j = 1:length(Mu)
    X_level_curve{j} = bvn_level_curve(Mu{j},Sigma,Level_curve,200);
    plot(X_level_curve{j}(:,1),X_level_curve{j}(:,2),'--',...
        'Color',ScatterColors{Class(j)},...
        'LineWidth',LineWidth)
    hold on
end
xlabel('X1')
ylabel('X2')
set(gca,'XLim',[-2 2],'YLim',[-2 2],'YTick',-2:2,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth,'Units',Units,'Position',[Axis_Left(5) Axis_Bottom(5) Axis_Width Axis_Height])

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
X_level_curve{j} = bvn_level_curve(Mu(j,1:2),Sigma(1:2,1:2,j),Level_curve,200);
end
obj = gmdistribution(Mu,Sigma,p);
[X,idx] = random(obj,n);
Y = Class(idx);
ax = subplot(2,3,6);
for c = 1:nclasses
    plot(X(Y==c,1),X(Y==c,2),'.',...
        'MarkerSize',MarkerSize,...
        'Color',ScatterColors{c})
    hold on
%plot(X(idx==j,1),X(Y==1,2),'.',X(Y==2,1),X(Y==2,2),'.',X(Y==3,1),X(Y==3,2),'.',...
%    X(Y==4,1),X(Y==4,2),'.',X_level_curve{1}(:,1),X_level_curve{1}(:,2),...
%    X_level_curve{2}(:,1),X_level_curve{2}(:,2),X_level_curve{3}(:,1),X_level_curve{3}(:,2),...
%    X_level_curve{4}(:,1),X_level_curve{4}(:,2))
end
for j = 1:J
    plot(X_level_curve{j}(:,1),X_level_curve{j}(:,2),'--',...
        'Color',ScatterColors{Class(j)},...
        'LineWidth',LineWidth)
    hold on
end
xlabel('X1')
ylabel('X2')
%ch = allchild(gca);
%for j = 1:J
%    set(ch(j),'MarkerSize',MarkerSize,'Color',Colors(Class(j)))
%end
%set(ch(6),'MarkerSize',MarkerSize,'Color','b')
%set(ch(7),'MarkerSize',MarkerSize,'Color','g')
%set(ch(8),'MarkerSize',MarkerSize,'Color','m')
set(gca,'XLim',[-4 4],'YLim',[-4 4],'YTick',-4:2:4,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth,'Units',Units,'Position',[Axis_Left(6) Axis_Bottom(6) Axis_Width Axis_Height])
hL = legend('Class 1','Class 2','Class 3','Class 4');
set(hL,'Units',Units,'Position',[Legend_Left(2) Legend_Bottom(2) Legend_Width(2) Legend_Height(2)],'Visible','On','Box',Box_Legend)
fname = '~/LOVEFest/Figures/Fig1_Lhat';
save_fig(gcf,fname)