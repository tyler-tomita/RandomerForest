clear
close all
clc

t = readtable('/cis/home/ttomita/Data/data_idx.txt','ReadVariableNames',false,'ReadRowNames',false);
dirnames = t.Var1;

Results = struct('Lhat',struct,'Time',struct);
Results2 = struct('Lhat',struct,'Time',struct);
Results3 = struct('Lhat',struct,'Time',struct);
j = 1;
for i = 1:length(dirnames)
    dirname = dirnames{i};
    fname = dir(strcat('/cis/home/ttomita/Data/',dirname,'/',dirname,'2*.mat'));
    fname2 = dir(strcat('/cis/home/ttomita/Data/',dirname,'/',dirname,'_Rotate.mat'));
    fname3 = dir(strcat('/cis/home/ttomita/Data/',dirname,'/',dirname,'_RerFdn.mat'));
        if isempty(fname) || isempty(fname2) || isempty(fname3)
        warning(sprintf('%s workspace file does not exist',dirname))
        rmidx(j) = i;
        j = j + 1;
    else
        load(fname.name,'Lhat','Time')
        Lhat.rerf = Lhat.f2;
        Lhat.rerfd = Lhat.f3;
        Lhat.rerfdr = Lhat.f4;
        Lhat = rmfield(Lhat,{'f2' 'f3' 'f4'});
        Time.rerf = Time.f2;
        Time.rerfd = Time.f3;
        Time.rerfdr = Time.f4;
        Time = rmfield(Time,{'f2' 'f3' 'f4'});
        Results.Lhat.(strrep(dirname,'-','__')) = Lhat;
        Results.Time.(strrep(dirname,'-','__')) = Time;
        load(fname2.name,'Lhat','Time')
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
        Results2.Lhat.(strrep(dirname,'-','__')) = Lhat;
        Results2.Time.(strrep(dirname,'-','__')) = Time;
        load(fname3.name,'Lhat','Time')
        Results3.Lhat.(strrep(dirname,'-','__')) = Lhat;
        Results3.Time.(strrep(dirname,'-','__')) = Time;
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
        Lhat_all.(clname)(i,:) = Results.Lhat.(strrep(dirname,'-','__')).(clname);
        Time_all.(clname)(i,:) = Results.Time.(strrep(dirname,'-','__')).(clname);
    end
    clnames2 = fieldnames(Results2.Lhat.(strrep(dirname,'-','__')));
    for j = 1:length(clnames2)
        clname = clnames2{j};
        Lhat_all.(clname)(i,:) = Results2.Lhat.(strrep(dirname,'-','__')).(clname)(end,:);
        Time_all.(clname)(i,:) = Results2.Time.(strrep(dirname,'-','__')).(clname);
    end
    clnames3 = fieldnames(Results3.Lhat.(strrep(dirname,'-','__')));
    for j = 1:length(clnames3)
        clname = clnames3{j};
        Lhat_all.(clname)(i,:) = Results3.Lhat.(strrep(dirname,'-','__')).(clname)(end,:);
        Time_all.(clname)(i,:) = Results3.Time.(strrep(dirname,'-','__')).(clname);
    end
end

clnames = cat(1,clnames,clnames2,clnames3);
ns = length(clnames);

for i = 1:ns
    clname = clnames{i};
    [Lhat_min.(clname),minidx.(clname)] = min(Lhat_all.(clname),[],2);
end

for i = 1:ns
    np(i) = length(Lhat_min.(clnames{i}));
end

if length(unique(np)) == 1
    np = unique(np);
    Lhat_ps = zeros(np,ns);
   
    for i = 1:ns
        Lhat_ps(:,i) = Lhat_min.(clnames{i});
    end
    
    Lhat_min_s = min(Lhat_ps,[],2);
    
    r_ps = zeros(size(Lhat_ps));
    
    for i = 1:ns
        r_ps(:,i) = Lhat_ps(:,i)./Lhat_min_s;
    end
    
    tau = 1:.1:ceil(max(max(r_ps)));
    
    rho_ps = zeros(length(tau),ns);
    for i = 1:length(tau)
        rho_ps(i,:) = sum(r_ps <= tau(i))/np;
    end
    
    Colors = linspecer(ns,'sequential');
    
    for i = 1:ns
        plot(tau,rho_ps(:,i),'Color',Colors(i,:))
        hold on
    end
    
    xlabel('r_p_,_s')
    ylabel('F(r_p_,_s)')
    legend(clnames)
end

fname = 'Performance_profile';
save_fig(gcf,fname)
    