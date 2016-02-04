clear
close all
clc

t = readtable('/cis/home/ttomita/Data/data_idx.txt','ReadVariableNames',false,'ReadRowNames',false);
dirnames = t.Var1;

Results = struct('Lhat',struct,'Time',struct);
j = 1;
for i = 1:length(dirnames)
    dirname = dirnames{i};
    fname = dir(strcat('/cis/home/ttomita/Data/',dirname,'/',dirname,'_final.mat'));
    fname2 = dir(strcat('/cis/home/ttomita/Data/',dirname,'/',dirname,'_rankRerF_final.mat'));
    if isempty(fname) || isempty(fname2)
        warning(sprintf('%s workspace file does not exist',dirname))
        j = j + 1;
    else
        load(fname2.name)
        Lhat2 = Lhat;
        Time2 = Time;
        load(fname.name)
        if length(fieldnames(Lhat)) == 4
            Lhat.rerfdnr = Lhat2.rerfdnr;
            Time.rerfdnr = Time2.rerfdnr;
            save(['/cis/home/ttomita/Data/',dirname,'/',strrep(fname.name,'.mat','2.mat')],'Lhat','Time','n','d','tr')
        else
            warning(sprintf('%s: not all classifiers evaluated',dirname))
        end
    end
end