clear
close all
clc

LineWidth = 2;
MarkerSize = 12;
FontSize = .2;
axWidth = 1.75;
axHeight = 1.75;
axLeft = [FontSize*4,FontSize*7+axWidth,FontSize*10+axWidth*2,...
    FontSize*13+axWidth*3];
axBottom = [FontSize*4,FontSize*4,FontSize*4,FontSize*4];
figWidth = axLeft(end) + axWidth + FontSize*4;
figHeight = axBottom(1) + axHeight + FontSize*4;
fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];

t = readtable('/cis/home/ttomita/Data/data_idx.txt','ReadVariableNames',false,'ReadRowNames',false);
dirnames = t.Var1;

Results = struct('Lhat',struct,'Time',struct);
j = 1;
for i = 1:length(dirnames)
    dirname = dirnames{i};
    fname = dir(strcat('/cis/home/ttomita/Data/',dirname,'/',dirname,'_2015_Nov_1.mat'));
    if isempty(fname)
        warning(sprintf('%s workspace file does not exist',dirname))
        rmidx(j) = i;
        j = j + 1;
    else
        load(fname.name,'Lhat','Time')
        if length(fieldnames(Lhat)) ~=4
            warning(sprintf('%s: not all classifiers evaluated',dirname))
            rmidx(j) = i;
            j = j + 1;
        else
            Results.Lhat.(strrep(dirname,'-','__')) = Lhat;
        end
    end
end

dirnames(rmidx) = [];

for i = 1:length(dirnames)
    dirname = dirnames{i};
    cd(dirname)
    if exist(strcat(dirname,'_train_R.dat'))
        FileName = strcat(dirname,'_train_R.dat');
    else
        FileName = strcat(dirname,'_R.dat');
    end
    X = dlmread(FileName,'\t',1,1);
    [n(i),d(i)] = size(X);
    cd ..
    clnames = fieldnames(Results.Lhat.(strrep(dirname,'-','__')));
    for j = 1:length(clnames)
        clname = clnames{j};
        Lhat_all.(clname)(i,:) = cat(2,mean(Results.Lhat.(strrep(dirname,'-','__')).(clname)(end,:,:),3),NaN(1,5-length(mean(Results.Lhat.(strrep(dirname,'-','__')).(clname)(end,:,:),3))));
    end
end

ns = length(clnames);

for i = 1:ns
    clname = clnames{i};
    [Lhat_min.(clname),minidx.(clname)] = min(Lhat_all.(clname),[],2);
end

for i = 1:ns
    np(i) = length(Lhat_min.(clnames{i}));
end

names = {'RF' 'RerF' 'RotRF' 'RerFd' 'RankRerFd'};
if length(unique(np)) == 1
    np = unique(np);
    Lhat_ps = zeros(np,ns);
   
    for i = 2:ns
        subplot(1,ns-1,i)
        plot(Lhat_min.(clnames{1}),Lhat_min.(clnames{i}),'.','MarkerSize',MarkerSize);
        hold on
        plot([0 1],[0 1],'-k','LineWidth',LineWidth)
        xlabel(sprintf('Error Rate (%s)',names{1}))
        ylabel(sprintf('Error Rate (%s)',names{i}))
        ax.LineWidth = LineWidth;
        ax.FontUnits = 'inches';
        ax.FontSize = FontSize;
        ax.Units = 'inches';
        ax.Position = [axLeft(i) axBottom(i) axWidth axHeight];
        ax.Box = 'off';
        ax.XLim = [0 1];
        ax.YLim = ax.XLim;
    end
end
    
save_fig(gcf,[rerfPath 'RandomerForest/Figures/Fig5_benchmark_scatter'])
