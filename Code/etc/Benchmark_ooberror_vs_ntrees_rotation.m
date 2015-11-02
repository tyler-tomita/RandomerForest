clear
close all
clc

t = readtable('/cis/home/ttomita/Data/data_idx.txt','ReadVariableNames',false,'ReadRowNames',false);
dirnames = t.Var1;

Results = struct('Lhat',struct,'Time',struct);
j = 1;
for i = 1:length(dirnames)
    dirname = dirnames{i};
    fname = dir(strcat('/cis/home/ttomita/Data/',dirname,'/',dirname,'_Rotate.mat'));
    if isempty(fname)
        warning(sprintf('%s workspace file does not exist',dirname))
        rmidx(j) = i;
        j = j + 1;
    else
        load(fname.name,'Lhat','Time')
        Lhat.rf_rot = Lhat.rf;
        Lhat.rerf_rot = Lhat.f2;
        Lhat.rerfd_rot = Lhat.f3;
        Lhat.rerfdn_rot = Lhat.f4;
        Lhat = rmfield(Lhat,{'rf' 'f2' 'f3' 'f4'});
        Time.rf_rot = Time.rf;
        Time.rerf_rot = Time.f2;
        Time.rerfd_rot = Time.f3;
        Time.rerfdn_rot = Time.f4;
        Time = rmfield(Time,{'rf' 'f2' 'f3' 'f4'});
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
    clnames = fieldnames(Results.Lhat.(strrep(dirname,'-','__')));
    for j = 1:length(clnames)
        clname = clnames{j};
        Lhat_all.(clname)(:,:,i) = Results.Lhat.(strrep(dirname,'-','__')).(clname);
        [~,minidx] = min(Lhat_all.(clname)(end,:,i),[],2);
        plot(Lhat_all.(clname)(:,minidx,i))
        hold on
    end
    xlabel('# trees')
    ylabel('oob error')
    title(dirname)
    legend(clnames)
    
    fname = strcat('/cis/home/ttomita/Data/',dirname,'/',dirname,'ooberror_vs_ntrees');
    save_fig(gcf,fname)
    
    close
end