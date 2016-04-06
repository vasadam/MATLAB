%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Read the statistics file, compute mean and MAE
rootdir = 'D:\NOAA\IGRA_PdT_from_2000';
files_total=0;
files_normal=0;
stationdirs = dir(rootdir);
for i=1:size(stationdirs,1)
    if (strcmp(stationdirs(i).name,'.') || strcmp(stationdirs(i).name,'..'))  % Skip '.' and '..'
        continue;
    end
    
    monthdirs = dir(fullfile(rootdir,stationdirs(i).name));    
    for j=1:size(monthdirs,1)        
        if (strcmp(monthdirs(j).name,'.') || strcmp(monthdirs(j).name,'..'))  % Skip '.' and '..'
            continue;
        end        
        
        hourfiles = dir(fullfile(rootdir,stationdirs(i).name,monthdirs(j).name));               
        for k=1:size(hourfiles)
            files_total = files_total+1;
            if (strcmp(hourfiles(k).name,'.') ...
                || strcmp(hourfiles(k).name,'..') ...
                || ~isempty(strfind(hourfiles(k).name,'stat')))  % Skip '.' and '..'
                continue;
            end
            
            fullfile(rootdir,stationdirs(i).name,monthdirs(j).name,hourfiles(k).name)
            STAT = load(fullfile(rootdir,stationdirs(i).name,monthdirs(j).name,hourfiles(k).name));
            avg = mean(STAT(:,2));
            stddev = std(STAT(:,2));
            if (size(STAT,1) > 3)
                normality = adtest(STAT(:,2));                
            else
                normality = 1;
            end
            files_normal = files_normal + (1-normality);
            
            filenameparts = strsplit(hourfiles(k).name,'.');
            fileoutname = strcat(filenameparts(1),'_stat.txt');
            fileout = char(fullfile(rootdir,stationdirs(i).name,monthdirs(j).name,fileoutname));
            dlmwrite(fileout, ...
                     [avg stddev normality], ...
                     'delimiter','\t');
        end
    end
end
