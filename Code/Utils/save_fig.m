function save_fig(h,fname)
    savefig(h,fname)
    saveas(h,fname,'pdf')
    saveas(h,fname,'png')
end