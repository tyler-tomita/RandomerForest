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
    catch
        warning(sprintf('%s workspace file does not exist',dirname))
        rmidx(j) = i;
	j = j + 1;
    end
end

dirnames(rmidx) = [];

for i = 1:length(dirnames)
    dirname = dirnames{i};
    clnames = fieldnames(Results.Lhat.(strrep(dirname,'-','__')));
    for j = 1:length(clnames)
        clname = clnames{j};
        Lhat_all.(clname)(i,:) = Results.Lhat.(strrep(dirname,'-','__')).(clname);
        Time_all.(clname)(i,:) = Results.Time.(strrep(dirname,'-','__')).(clname);
    end
end

plot(Lhat_all.rf(:,1),Lhat_all.f4(:,1),'b.');
hold on
plot([0 1],[0 1],'-k')

lspec = {'bs','rs','gs' 'cs' 'ms' 'ks'};
facespec = {'b','r','g' 'c' 'm' 'k'};

save_fig(gcf,'/cis/home/ttomita/Data/Fig4_Real_Data_panel_b')
