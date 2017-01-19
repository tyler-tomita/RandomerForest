function save_fig(h,fname,FileType)
    if iscell(FileType)
        FileType = unique(FileType);
        FileType = FileType(ismember(FileType,{'fig','pdf','png'}));
        for i = 1:length(FileType)
            if strcmp(FileType{i},'fig')
                savefig(h,fname)
            else
                saveas(h,fname,FileType{i})
            end
        end
    else
        if ismember(FileType,{'fig','pdf','png'})
            if strcmp(FileType,'fig')
                savefig(h,fname)
            else
                saveas(h,fname,FileType)
            end
        end
    end
end