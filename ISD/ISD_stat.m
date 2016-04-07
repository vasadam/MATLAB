newline = java.lang.System.getProperty('line.separator');
% hourFile = fullfile('D:\NOAA\IGRA_PerStationStatistics',stationdirs(i).name,hourFiles(j).name);
ISDStatisticsRootDir = 'D:\NOAA\ISD_parsed';
IGRAStatisticsRootDir = 'D:\NOAA\IGRA_PerStationStatistics';
ISDStatisticsDestDir = 'D:\NOAA\ISD_stat';
stationFiles = dir(ISDStatisticsRootDir);
% for n=1:size(stationFiles,1)
parfor (n=1:size(stationFiles,1),4)
    if (strcmp(stationFiles(n).name,'.') || strcmp(stationFiles(n).name,'..'))  % Skip '.' and '..'
        continue;
    end         
    
    fileNameParts = strsplit(stationFiles(n).name,'.'); 
    stationName = fileNameParts{1};    
    stationFile = fopen(fullfile(ISDStatisticsRootDir,strcat(stationName,'.txt')),'r');
    line = fgetl(stationFile);
    measurements = cell(1,3);
    k = 0;
    while(ischar(line))
        k = k+1;
        values = textscan(line,'%s %s %s %f %f\n');
        measurements(k,:) = {datetime(strcat(values{2},values{3}),'InputFormat','yyyyMMddHHmm'),values{4},values{5}};    
        line = fgetl(stationFile);
    end
    fclose(stationFile);

    % Hour differences to take into account during the derivation
    hourDiffsMatter = [1 3 6 9 12 24];           

    % Indices of columns in the measurements matrix containing the derivatives
    pDerColIndices = cell(1,4);
    teDerColIndices = cell(1,4);               
    derColIndicesStartingValues = size(measurements,2)+1 : size(measurements,2) + size(hourDiffsMatter,2); % Starting values of column indices of derivatives             
    for i=1:4
        pDerColIndices{i} = containers.Map(hourDiffsMatter, derColIndicesStartingValues + (2*i-2)*size(hourDiffsMatter,2));                      
        teDerColIndices{i} = containers.Map(hourDiffsMatter, derColIndicesStartingValues + (2*i-1)*size(hourDiffsMatter,2));    
    end

    %% 1st order derivatives                                    
    for i=2:size(measurements,1)-1    
        neighborIndices = containers.Map('KeyType','int32','ValueType','uint64');
        pDerivatives1st = containers.Map('KeyType','uint32','ValueType','double');
        teDerivatives1st = containers.Map('KeyType','uint32','ValueType','double');

        % Find significant measurements backward
        for j=i-1:-1:1
            hourDiff = hours(measurements{j,1} - measurements{i,1});
            if (-hourDiff > max(hourDiffsMatter))
                break;
            end
            if (~ismember(-hourDiff,hourDiffsMatter))
                continue;
            end
            neighborIndices(hourDiff) = j;
        end
        % Find significant measurements forward
        for j=i+1:size(measurements,1)
            hourDiff = hours(measurements{j,1} - measurements{i,1});
            if (hourDiff > max(hourDiffsMatter))
                break;
            end
            if (~ismember(hourDiff,hourDiffsMatter))
                continue;
            end
            neighborIndices(hourDiff) = j;
        end

        % Apply central differences method
        if (isempty(neighborIndices))
            continue;
        end
        keySet = cell2mat(keys(neighborIndices));
        for j=size(keySet,2):-1:1
            % It is sufficient to iterate only over positive key values.        
            hourDiff = keySet(j);
            if (hourDiff<0)
                break;
            end
            % Skip keys which doesn't have a negative pair.        
            if (~ismember(-hourDiff,keySet))
                continue;
            end
            % Apply central differences method for P
            pNext = measurements{neighborIndices(hourDiff),2};
            pPrev = measurements{neighborIndices(-hourDiff),2};
            pDerivatives1st(hourDiff) = (pNext-pPrev) / (2*double(hourDiff));

            % Apply central differences method for Te
            tePrev = measurements{neighborIndices(-hourDiff),3};
            teNext = measurements{neighborIndices(hourDiff),3};        
            % Te value may be missing, in that case skip calculation
            if (teNext == -99999 || tePrev == -99999)
                continue;
            end
            teDerivatives1st(hourDiff) = (teNext-tePrev) / (2*double(hourDiff));                            
        end

        % Write results into the matrix
        keySet = cell2mat(keys(pDerivatives1st));
        colIndices = pDerColIndices{1}; 
        for j=1:size(keySet,2)        
            measurements{i,colIndices(keySet(j))} = pDerivatives1st(keySet(j));
        end    

        keySet = cell2mat(keys(teDerivatives1st));
        colIndices = teDerColIndices{1};    
        for j=1:size(keySet,2)
            measurements{i,colIndices(keySet(j))} = teDerivatives1st(keySet(j));
        end
    end

    %% 2nd order derivatives 
    for i=3:size(measurements,1)-2    
        for j=1:size(hourDiffsMatter,2)
            hourDiff = hourDiffsMatter(j);
            pPrev = measurements{i-1,pDerColIndices{1}(hourDiff)};
            pNext = measurements{i+1,pDerColIndices{1}(hourDiff)};
            tePrev = measurements{i-1,teDerColIndices{1}(hourDiff)};
            teNext = measurements{i+1,teDerColIndices{1}(hourDiff)};
            if (~isempty(pPrev) && ~isempty(pNext))        
                measurements{i,pDerColIndices{2}(hourDiff)} = (pNext-pPrev) / (2*double(hourDiff));
            end
            if (~isempty(tePrev) && ~isempty(teNext))                
                measurements{i,teDerColIndices{2}(hourDiff)} = (teNext-tePrev) / (2*double(hourDiff));
            end        
        end
    end

    %% 3rd order derivatives 
    for i=4:size(measurements,1)-3    
        for j=1:size(hourDiffsMatter,2)
            hourDiff = hourDiffsMatter(j);
            pPrev = measurements{i-1,pDerColIndices{2}(hourDiff)};
            pNext = measurements{i+1,pDerColIndices{2}(hourDiff)};
            tePrev = measurements{i-1,teDerColIndices{2}(hourDiff)};
            teNext = measurements{i+1,teDerColIndices{2}(hourDiff)};
            if (~isempty(pPrev) && ~isempty(pNext))        
                measurements{i,pDerColIndices{3}(hourDiff)} = (pNext-pPrev) / (2*double(hourDiff));
            end
            if (~isempty(tePrev) && ~isempty(teNext))                
                measurements{i,teDerColIndices{3}(hourDiff)} = (teNext-tePrev) / (2*double(hourDiff));
            end        
        end
    end

    %% 4th order derivatives 
    for i=5:size(measurements,1)-4    
        for j=1:size(hourDiffsMatter,2)
            hourDiff = hourDiffsMatter(j);
            pPrev = measurements{i-1,pDerColIndices{3}(hourDiff)};
            pNext = measurements{i+1,pDerColIndices{3}(hourDiff)};
            tePrev = measurements{i-1,teDerColIndices{3}(hourDiff)};
            teNext = measurements{i+1,teDerColIndices{3}(hourDiff)};
            if (~isempty(pPrev) && ~isempty(pNext))        
                measurements{i,pDerColIndices{4}(hourDiff)} = (pNext-pPrev) / (2*double(hourDiff));
            end
            if (~isempty(tePrev) && ~isempty(teNext))                
                measurements{i,teDerColIndices{4}(hourDiff)} = (teNext-tePrev) / (2*double(hourDiff));
            end        
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Read IGRA measurements and insert derivative values from ISD    
    IGRAHours = {'00','12'};
    for hour = IGRAHours
        hourFileName = fullfile(IGRAStatisticsRootDir,stationName,strcat(hour{1},'.txt'));
        if (exist(hourFileName,'file') == 0)
            continue;
        end
        hourFile = fopen(hourFileName,'r');
        values_new = cell(1);
        line = fgetl(hourFile);
        k = 0;
        j_prev = 1;
        while(ischar(line))    
            values = textscan(line,'%s %s %s %s %s %s %s %s %s %s %s %s\n');
            year = values{1}{1};
            month = values{2}{1};
            if (size(month)==1)
                month = strcat('0',month);
            end
            day = values{3}{1};
            if (size(day)==1)
                day = strcat('0',day);
            end            
%             fileNameParts = strsplit(hourFiles(j).name,'.');
%             hour = fileNameParts{1};    
            minute = '00';  % Its value is always 00
        %     zLL = values{4};
        %     pLL = values{5};
        %     tLL = values{6};
        %     rhLL = values{7};
        %     tvLL = values{8};
        %     teLL = values{9};
        %     z500 = values{10};   
        %     z500_calc = values{11};            
        %     E = values{12};

            % Find the corresponding measurement in the ISD measurements cell array
            dt = datetime(strcat(year,month,day,hour{1},minute),'InputFormat','yyyyMMddHHmm');
            for j=j_prev:size(measurements,1)
                hourDiff = hours(dt - measurements{j,1});
                if (hourDiff < 0) % Stepped over the current datetime -> no match
                    j_prev = j;
                    break;
                end
                if (hourDiff == 0) % Found match
                    k = k+1;
                    values_new{k} = cellstr([strrep(line,' ',','),',',...
                                             strjoin(cellfun(@(x) num2str(x),...
                                                             measurements(j,derColIndicesStartingValues(1):size(measurements,2)),...
                                                             'UniformOutput',false),...
                                                     ',')]); 
                    j_prev = j;
                    break;
                end
            end    
            line = fgetl(hourFile);
        end
        fclose(hourFile);

        %% TODO
        mkdir(fullfile(ISDStatisticsDestDir,stationName));
        stationFile = fopen(fullfile(ISDStatisticsDestDir,stationName,strcat(hour{1},'.txt')),'w');
        for k=1:size(values_new,2)
            k
            fprintf(stationFile,'%s',[char(values_new{k}),char(newline)]);
        end
        fclose(stationFile);   
    end        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Import data from the recently created text files
%  and perform the statistical analysis on them
stationDirs = dir(ISDStatisticsDestDir);
for n=1:size(stationDirs,1)
    if (strcmp(stationDirs(n).name,'.') || strcmp(stationDirs(n).name,'..'))  % Skip '.' and '..'
        continue;
    end   
    
    hourFiles = dir(fullfile(ISDStatisticsDestDir,stationDirs(n).name));
    for m=1:size(hourFiles,1)
        if (strcmp(hourFiles(n).name,'.') || strcmp(hourFiles(n).name,'..'))  % Skip '.' and '..'
            continue;
        end                
        
        formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
        hourFileID = fopen(hourFiles(n).name,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', ',', 'EmptyValue', NaN, 'ReturnOnError', false);
        fclose(hourFileID);
        
        % Convert cell array to matrix
        dataArray = [dataArray{1:end-1}];
    end
end