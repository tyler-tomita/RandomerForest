clear
close all
clc

t = readtable('/cis/home/ttomita/Data/sorted_list.txt','ReadVariableNames',false,'ReadRowNames',false);
dirnames = t.Var1;

Results = struct('Lhat',struct,'Time',struct,'sp',struct);
Results2 = struct('Lhat',struct,'Time',struct,'sp',struct);
j = 1;
for i = 1:length(dirnames)
    dirname = dirnames{i};
    fname = dir(strcat('/cis/home/ttomita/Data/',dirname,'/',dirname,'_RerFr.mat'));
    fname2 = dir(strcat('/cis/home/ttomita/Data/',dirname,'/',dirname,'_RF.mat'));
    if isempty(fname) 
        fprintf('%s\n%s\n',fname.name)
        warning(sprintf('%s workspace file does not exist',dirname))
        rmidx(j) = i;
	j = j + 1;
    elseif isempty(fname2)
        fprintf('%s\n%s\n',fname2.name)
        warning(sprintf('%s workspace file does not exist',dirname))
        rmidx(j) = i;
        j = j + 1;
    else
        load(fname.name,'Lhat','Time','sp')
        Results.Lhat.(strrep(dirname,'-','__')) = Lhat;
        Results.Time.(strrep(dirname,'-','__')) = Time;
	Results.sp.(strrep(dirname,'-','__')) = sp;
        load(fname2.name,'Lhat','Time','sp')
        Results2.Lhat.(strrep(dirname,'-','__')) = Lhat;
        Results2.Time.(strrep(dirname,'-','__')) = Time;
	Results2.sp.(strrep(dirname,'-','__')) = sp;
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
    Lhat_all.f5(:,:,i) = Results.Lhat.(strrep(dirname,'-','__')).f5;
    Time_all.f5(i,:) = Results.Time.(strrep(dirname,'-','__')).f5;
    sp_all.f5(i,:) = Results.sp.(strrep(dirname,'-','__')).f5;
    Lhat_all.rf(:,:,i) = Results2.Lhat.(strrep(dirname,'-','__')).rf;
    Time_all.rf(i,:) = Results2.Time.(strrep(dirname,'-','__')).rf;
    sp_all.rf(i,:) = Results2.sp.(strrep(dirname,'-','__')).rf;
end

lspec = {'bs','rs','gs' 'cs' 'ms'};
facespec = {'b','r','g' 'c' 'm'};


%%%%%change so that take lowest Lhat of the three and sum time over all ks%%%%%%%%%%%%
clnames = {'f5','rf'};
for i = 1:length(clnames)
    clname = clnames{i};
    [Lhat_min.(clname),minidx.(clname)] = min(Lhat_all.(clname)(end,:,:),[],2);
    Lhat_mean.(clname) = mean(Lhat_min.(clname));
    Lhat_sem.(clname) = std(Lhat_min.(clname))/sqrt(length(dirnames));
    Time_mean.(clname) = mean(mean(Time_all.(clname),2));
    Time_sem.(clname) = std(mean(Time_all.(clname),2))/sqrt(length(dirnames));
    for j = 1:length(dirnames)
        sp_min.(clname)(j,1) = sp_all.(clname)(j,minidx.(clname)(j));
    end
end

scatter(sp_min.rf,sp_min.f5);
mxrf = max(sp_min.rf);
mxf5 = max(sp_min.f5);
mx = max([mxrf;mxf5]);
hold on
plot([0 mx],[0 mx],'-k')
xlabel('RF')
ylabel('RerF(r)')
title('Complexity')
save_fig(gcf,'/cis/home/ttomita/Data/Fig4_Complexity')
