newline = java.lang.System.getProperty('line.separator');
months_of_seasons = containers.Map;
months_of_seasons('winter') = [12 1 2];
months_of_seasons('spring') = [3 4 5];
months_of_seasons('summer') = [6 7 8];
months_of_seasons('autumn') = [9 10 11];
ISDStatisticsRootDir = 'D:\NOAA\ISD_stat_tend';
ISDStatisticsDestDir = 'D:\NOAA\ISD_stat_tend_corr_per_season';
stationDirs = dir(ISDStatisticsRootDir);
parfor i=1:size(stationDirs,1)
    if (strcmp(stationDirs(i).name,'.') || strcmp(stationDirs(i).name,'..'))  % Skip '.' and '..'
        continue;
    end     
    
    hourFiles = dir(fullfile(ISDStatisticsRootDir,stationDirs(i).name));
    for j=1:size(hourFiles,1)
        if (strcmp(hourFiles(j).name,'.') || strcmp(hourFiles(j).name,'..'))  % Skip '.' and '..'
            continue;
        end      
        
        %% Format string for each line of text:
        % For more information, see the TEXTSCAN documentation.
        formatSpec = '%f%f%f%f%f%*s%*s%*s%f%*s%*s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';        
        %% Open the text file.
        fileID = fopen(fullfile(ISDStatisticsRootDir,stationDirs(i).name,hourFiles(j).name),'r');

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
        elevation = dataArray{:, 4};
        P = dataArray{:, 5};
        Te = dataArray{:, 6};
        E = dataArray{:, 7};
        d1P1h = dataArray{:, 8};
        d1P3h = dataArray{:, 9};
        d1P6h = dataArray{:, 10};
        d1P9h = dataArray{:, 11};
        d1P12h = dataArray{:, 12};
        d1P24h = dataArray{:, 13};
        d1Te1h = dataArray{:, 14};
        d1Te3h = dataArray{:, 15};
        d1Te6h = dataArray{:, 16};
        d1Te9h = dataArray{:, 17};
        d1Te12h = dataArray{:, 18};
        d1Te24h = dataArray{:, 19};
        d2P1h = dataArray{:, 20};
        d2P3h = dataArray{:, 21};
        d2P6h = dataArray{:, 22};
        d2P9h = dataArray{:, 23};
        d2P12h = dataArray{:, 24};
        d2P24h = dataArray{:, 25};
        d2Te1h = dataArray{:, 26};
        d2Te3h = dataArray{:, 27};
        d2Te6h = dataArray{:, 28};
        d2Te9h = dataArray{:, 29};
        d2Te12h = dataArray{:, 30};
        d2Te24h = dataArray{:, 31};
        d3P1h = dataArray{:, 32};
        d3P3h = dataArray{:, 33};
        d3P6h = dataArray{:, 34};
        d3P9h = dataArray{:, 35};
        d3P12h = dataArray{:, 36};
        d3P24h = dataArray{:, 37};
        d3Te1h = dataArray{:, 38};
        d3Te3h = dataArray{:, 39};
        d3Te6h = dataArray{:, 40};
        d3Te9h = dataArray{:, 41};
        d3Te12h = dataArray{:, 42};
        d3Te24h = dataArray{:, 43};
        d4P1h = dataArray{:, 44};
        d4P3h = dataArray{:, 45};
        d4P6h = dataArray{:, 46};
        d4P9h = dataArray{:, 47};
        d4P12h = dataArray{:, 48};
        d4P24h = dataArray{:, 49};
        d4Te1h = dataArray{:, 50};
        d4Te3h = dataArray{:, 51};
        d4Te6h = dataArray{:, 52};
        d4Te9h = dataArray{:, 53};
        d4Te12h = dataArray{:, 54};
        d4Te24h = dataArray{:, 55};         
        
        corrs = containers.Map('KeyType','char','ValueType','double');
        pVals = containers.Map('KeyType','char','ValueType','double');        
        %% Calculate correlation coefficients    
        for k=1:size(dataArray,2)
            k
            % Skip year, month, day, elevation and E
            if (k==1 || k==2 || k==3 || k==4 || k==7)
                continue;
            end
            values = dataArray{:,k};
            if (~isempty(values(~isnan(values))))
                %% same-day correlation coefficient
                for season = keys(months_of_seasons)                
                    condition = ~isnan(values) & ismember(month,months_of_seasons(char(season))); 
                    key = [char(season),' ',num2str(k),' ','0'];
                    if (~isempty(values(condition)))                        
                        [corrs(key), pVals(key)] = corr(E(condition), values(condition));
                    else
                        corrs(key) = NaN;
                        pVals(key) = NaN;
                    end
                end    
                %% previous-days correlation coefficients
                for d=-1:-1:-5
                    month_tmp = [];
                    prev_values_tmp = [];
                    E_tmp = [];
                    n=0;   
                    for q=(1-d):size(E)
                        % check if previous day measurement is available                        
                        fileNameParts = strsplit(hourFiles(j).name,'.');
                        hour = fileNameParts{1};    
                        minute = '00';  % Its value is always 00
                        dtNow = datetime(strcat(num2str(year(q)),...
                                                num2str(month(q),'%02d'),...
                                                num2str(day(q),'%02d'),...
                                                num2str(str2double(hour),'%02d'),...
                                                minute),...
                                         'InputFormat', 'yyyyMMddHHmm');
                        for r=q-1:-1:1
                            dtPrev = datetime(strcat(num2str(year(r)),...
                                                     num2str(month(r),'%02d'),...
                                                     num2str(day(r),'%02d'),...
                                                     num2str(str2double(hour),'%02d'),...
                                                     minute),...
                                              'InputFormat', 'yyyyMMddHHmm');
                            hourDiff = hours(dtPrev - dtNow);
                            if (hourDiff < d*24)
                                break;
                            end
                            % If available (and not NaN), store E and
                            % measurement in another array
                            if (hourDiff == d*24)
                                if (isnan(values(r)))
                                    break;
                                end
                                n=n+1;
                                month_tmp(n) = month(q);
                                E_tmp(n) = E(q);
                                prev_values_tmp(n) = values(r);
                            end
                        end                        
                    end   
                    
                    % per-season previous-day correlation coefficients
                    for season = keys(months_of_seasons)
                        condition = ismember(month_tmp,months_of_seasons(char(season)));  
                        key = [char(season),' ',num2str(k),' ',num2str(d)];
                        if (~isempty(E_tmp(condition)))                            
                            [corrs(key), pVals(key)] = corr(E_tmp(condition)', prev_values_tmp(condition)');
                        else
                            corrs(key) = NaN;
                            pVals(key) = NaN;                            
                        end
                    end                    
                end                
            else
                % If column is empty, set all coefficients to NaN
                for d=0:-1:-5
                    for season = keys(months_of_seasons)
                        key = [char(season),' ',num2str(k),' ',num2str(d)];
                        corrs(key) = NaN;
                        pVals(key) = NaN;                        
                    end
                end                
            end
        end
        
        %% Write corr,pval values to file
        mkdir(fullfile(ISDStatisticsDestDir,stationDirs(i).name));
        hourFile = fopen(fullfile(ISDStatisticsDestDir,stationDirs(i).name,hourFiles(j).name),'w');
        for season = keys(months_of_seasons)
            for k=1:size(dataArray,2)
                % Skip year, month, day, elevation and E
                if (k==1 || k==2 || k==3 || k==4 || k==7)
                    continue;
                end
                for d=0:-1:-5
                    fprintf(hourFile,'%s',[char(season),...
                                           ',',...
                                           num2str(k),...
                                           ',',...
                                           num2str(d),...
                                           ',',...                                       
                                           num2str(corrs([char(season),' ',num2str(k),' ',num2str(d)])),...
                                           ',',...
                                           num2str(pVals([char(season),' ',num2str(k),' ',num2str(d)])),...
                                           char(newline)]);
                end 
            end
        end
        fclose(hourFile); 
    end
end
