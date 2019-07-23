% a. clean
if exist('../output/ConSav','dir')
   rmdir('../output/ConSav','s'); 
end
mkdir('../output/ConSav');

% b. figures
dirs = dir('figs_tabs\*');
for i = 3:numel(dirs)
    
    folder = dirs(i).name;
    files = dir(sprintf('figs_tabs/%s/*.pdf',folder));
    
    for j = 1:numel(files)
        fn = sprintf('figs_tabs/%s/%s',folder,files(j).name);
        fnout = sprintf('../output/ConSav/%s_%s',folder,files(j).name);
        copyfile(fn,fnout);
    end
   
end

% c. tables
copyfile('figs_tabs/main_table.tex','../output/ConSav/main_table.tex');