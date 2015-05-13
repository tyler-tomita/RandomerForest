clear
close all
clc

t = readtable('/cis/home/ttomita/Data/data_idx.txt','ReadVariableNames',false,'ReadRowNames',false);
dirnames = t.Var1;
Results = struct('Lhat',struct,'Time',struct);
j = 1;
for i = 1:length(dirnames)
    dirname = dirnames{i};
    try
        load(strcat('/cis/home/ttomita/Data/',dirname,'/',dirname,'.mat'),'Lhat','Time')
        Results.Lhat.(strrep(dirname,'-','__')) = Lhat;
        Results.Time.(strrep(dirname,'-','__')) = Time;
	load(strcat('/cis/home/ttomita/Data/',dirname,'/',dirname,'_NB.mat'),'Lhat')
    catch
        warning(sprintf('%s workspace file does not exist',dirname))
        rmidx(j) = i;
	j = j + 1;
    end
end

dirnames(rmidx) = [];
Lhat_nb = NaN(length(dirnames),1);

for i = 1:length(dirnames)
    dirname = dirnames{i};
    load(strcat('/cis/home/ttomita/Data/',dirname,'/',dirname,'_NB.mat'),'Lhat')
    Lhat_nb(i) = Lhat;
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
        Lhat_all.(clname)(i,:) = Results.Lhat.(strrep(dirname,'-','__')).(clname);
        Time_all.(clname)(i,:) = Results.Time.(strrep(dirname,'-','__')).(clname);
    end
end

lspec = {'bs','rs','gs' 'cs' 'ms' 'ks'};
facespec = {'b','r','g' 'c' 'm' 'k'};

Lhat_delta = Lhat_all.f4(:,end) - Lhat_all.rf(:,1);

%subplot(1,4,1)
figure(1)
plot(d,Lhat_delta,'.');
xlabel('D')
ylabel('Lhat(R''erF(s+d+r)) - Lhat(RF)')
save_fig(gcf,'/cis/home/ttomita/Data/Fig5_panel_1')

%subplot(1,4,2)
figure(2)
plot(n,Lhat_delta,'.')
xlabel('n')
ylabel('Lhat(R''erF(s+d+r)) - Lhat(RF)')
save_fig(gcf,'/cis/home/ttomita/Data/Fig5_panel_2')

%subplot(1,4,3)
figure(3)
plot(n./d,Lhat_delta,'.')
xlabel('n/d')
ylabel('Lhat(R''erF(s+d+r)) - Lhat(RF)')
save_fig(gcf,'/cis/home/ttomita/Data/Fig5_panel_3')

%subplot(1,4,4)
figure(4)
plot(Lhat_nb,Lhat_delta,'.')
xlabel('Lhat(Naive Bayes)')
ylabel('Lhat(R''erF(s+d+r)) - Lhat(RF)')
save_fig(gcf,'/cis/home/ttomita/Data/Fig5_panel_4')
