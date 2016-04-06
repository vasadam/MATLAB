rootDir = 'D:\NOAA\ISD_stat';
destDir = 'D:\NOAA\ISD_stat_bin';
stationDirs = dir(rootDir);
parfor i=1:size(stationDirs,1)
    if (strcmp(stationDirs(i).name,'.') || strcmp(stationDirs(i).name,'..'))  % Skip '.' and '..'
        continue;
    end         
    
    hourFiles = dir(fullfile(rootDir,stationDirs(i).name));
    for j=1:size(hourFiles,1)
        if (strcmp(hourFiles(j).name,'.')...
            || strcmp(hourFiles(j).name,'..')...
            || isdir(hourFiles(j).name))  % Skip '.', '..' and folders
            continue;
        end      
        
        %% Format string for each line of text:
        % For more information, see the TEXTSCAN documentation.
        formatSpec = '%f%f%f%f%f%*s%*s%*s%f%*s%*s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';        
        %% Open the text file.
        fileID = fopen(fullfile(rootDir,stationDirs(i).name,hourFiles(j).name),'r');

        %% Read columns of data according to format string.
        % This call is based on the structure of the file used to generate this
        % code. If an error occurs for a different file, try regenerating the code
        % from the Import Tool.
        dataArray = textscan(fileID, formatSpec, 'Delimiter', ',', 'EmptyValue', NaN, 'ReturnOnError', false);
        dataArray = dataArray(1:end-1); % delete the last, empty column

        %% Close the text file.
        fclose(fileID);

        %% Allocate imported array to column variable names
        year = dataArray{:, 1};
        month = dataArray{:, 2};
        day = dataArray{:, 3};
        fileNameParts = strsplit(hourFiles(j).name,'.');
        hour = fileNameParts{1};            
        
        %% Convert columns 1-3 to datetime
        dt = datetime.empty;
        for q=1:size(dataArray{1},1)
            dt(q) = datetime(strcat(num2str(year(q)),...
                                    num2str(month(q),'%02d'),...
                                    num2str(day(q),'%02d'),...
                                    num2str(str2double(hour),'%02d'),...
                                    '00'),...
                             'InputFormat', 'yyyyMMddHHmm');     
        end
        dataArray{3} = dt';
        dataArray = dataArray(3:end);
        
        %% Write variables to binary files
        destStationDirName = fullfile(destDir,stationDirs(i).name);
        if (~isdir(destStationDirName))
            mkdir(destStationDirName);
        end
        parsave(fullfile(destStationDirName,[hour '.bin']), dataArray);
    end
end
