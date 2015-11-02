clear
close all
clc

t = readtable('/cis/home/ttomita/Data/sorted_list.txt','ReadVariableNames',false,'ReadRowNames',false);
dirnames = t.Var1;

Results = struct('Lhat',struct,'Time',struct,'sp',struct);
j = 1;
for i = 1:length(dirnames)
    dirname = dirnames{i};
    fname = dir(strcat('/cis/home/ttomita/Data/',dirname,'/',dirname,'_RerFr.mat'));
    if isempty(fname)
        fprintf('%s\n%s\n',fname.name)
        warning(sprintf('%s workspace file does not exist',dirname))
        rmidx(j) = i;
	j = j + 1;
    else
        load(fname.name,'Lhat','Time','sp')
        Results.Lhat.(strrep(dirname,'-','__')) = Lhat;
        Results.Time.(strrep(dirname,'-','__')) = Time;
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
    clnames = {'f5'};
    for j = 1:length(clnames)
        clname = clnames{j};
        Lhat_all.(clname)(:,:,i) = Results.Lhat.(strrep(dirname,'-','__')).(clname);
        Time_all.(clname)(i,:) = Results.Time.(strrep(dirname,'-','__')).(clname);
    end
end

lspec = {'bs','rs','gs' 'cs' 'ms'};
facespec = {'b','r','g' 'c' 'm'};


%%%%%change so that take lowest Lhat of the three and sum time over all ks%%%%%%%%%%%%
clnames = {'f5'};
for i = 1:length(clnames)
    clname = clnames{i};
    [Lhat_min.(clname),minidx.(clname)] = min(Lhat_all.(clname)(end,:,:),[],2);
    Lhat_mean.(clname) = mean(Lhat_min.(clname));
    Lhat_sem.(clname) = std(Lhat_min.(clname))/sqrt(length(dirnames));
    Time_mean.(clname) = mean(mean(Time_all.(clname),2));
    Time_sem.(clname) = std(mean(Time_all.(clname),2))/sqrt(length(dirnames));
    plot(Time_mean.(clname),Lhat_mean.(clname),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
    hold on
end

xlabel('Training Time (sec)')
ylabel('Lhat')
legend('RerF(r)')

for i = 1:length(clnames)
    clname = clnames{i};
    Mu = [Time_mean.(clname) Lhat_mean.(clname)];
    Sigma = cov(mean(Time_all.(clname),2),squeeze(Lhat_min.(clname)));
    X_level_curve = bvn_level_curve(Mu,Sigma,0.1,200);
    plot(X_level_curve(:,1),X_level_curve(:,2),'--',...
         'Color',facespec{i},...
         'LineWidth',1.5)
    hold on
end

save('/cis/home/ttomita/Data/Fig4_RerFr.mat','Results','Lhat_mean','Lhat_sem','Time_mean','Time_sem')
save_fig(gcf,'/cis/home/ttomita/Data/Fig4_Real_Data_Panel_A_RerFr')
