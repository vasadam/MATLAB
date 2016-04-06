newline = java.lang.System.getProperty('line.separator');
% hourFile = fullfile('D:\NOAA\IGRA_PerStationStatistics',stationdirs(i).name,hourFiles(j).name);
ISDStatisticsRootDir = 'D:\NOAA\ISD_parsed';
IGRAStatisticsRootDir = 'D:\NOAA\IGRA_PerStationStatistics';
ISDStatisticsDestDir = 'D:\NOAA\ISD_stat_tend';
stationFiles = dir(ISDStatisticsRootDir);
% for n=1:size(stationFiles,1)
for n=1:size(stationFiles,1)
    if (strcmp(stationFiles(n).name,'.') || strcmp(stationFiles(n).name,'..'))  % Skip '.' and '..'
        continue;
    end         
    stationFiles(n).name
    
    if (~strcmp(stationFiles(n).name,'NLM00006235.txt'))
        continue;
    end
    
    fileNameParts = strsplit(stationFiles(n).name,'.'); 
    stationName = fileNameParts{1};    
    stationFile = fopen(fullfile(ISDStatisticsRootDir,strcat(stationName,'.txt')),'r');
    line = fgetl(stationFile);
    measurements = cell(1,3);
    k = 0;
    while(ischar(line))
        k = k+1
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
        i
        neighborIndices = containers.Map('KeyType','int32','ValueType','uint64');
        pDerivatives1st = containers.Map('KeyType','uint32','ValueType','double');
        teDerivatives1st = containers.Map('KeyType','uint32','ValueType','double');

        % Find significant measurements backward
        for j=i-1:-1:1
            hourDiff = hours(measurements{i,1} - measurements{j,1});
            if (hourDiff > max(hourDiffsMatter))
                break;
            end
            if (~ismember(hourDiff,hourDiffsMatter))
                continue;
            end
            neighborIndices(hourDiff) = j;
        end

        % Apply backward differences method
        if (isempty(neighborIndices))
            continue;
        end        
        keySet = cell2mat(keys(neighborIndices));
        for j=1:size(keySet,2)
            hourDiff = keySet(j);
            % Apply backward differences method for P
            pPrev = measurements{neighborIndices(hourDiff),2};
            pNow = measurements{i,2};
            pDerivatives1st(hourDiff) = (pNow-pPrev) / double(hourDiff);

            % Apply backward differences method for Te
            tePrev = measurements{neighborIndices(hourDiff),3};     
            teNow = measurements{i,3};            
            % Te value may be missing, in that case skip calculation
            if (teNow == -99999 || tePrev == -99999)
                continue;
            end
            teDerivatives1st(hourDiff) = (teNow-tePrev) / double(hourDiff);                            
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
        i
        for j=1:size(hourDiffsMatter,2)
            hourDiff = hourDiffsMatter(j);
            pPrev = measurements{i-1,pDerColIndices{1}(hourDiff)};
            pNow = measurements{i,pDerColIndices{1}(hourDiff)};
            tePrev = measurements{i-1,teDerColIndices{1}(hourDiff)};
            teNow = measurements{i,teDerColIndices{1}(hourDiff)};
            if (~isempty(pPrev) && ~isempty(pNow))        
                measurements{i,pDerColIndices{2}(hourDiff)} = (pNow-pPrev) / double(hourDiff);
            end
            if (~isempty(tePrev) && ~isempty(teNow))                
                measurements{i,teDerColIndices{2}(hourDiff)} = (teNow-tePrev) / double(hourDiff);
            end        
        end
    end

    %% 3rd order derivatives 
    for i=4:size(measurements,1)-3    
        i
        for j=1:size(hourDiffsMatter,2)
            hourDiff = hourDiffsMatter(j);
            pPrev = measurements{i-1,pDerColIndices{2}(hourDiff)};
            pNow = measurements{i,pDerColIndices{2}(hourDiff)};
            tePrev = measurements{i-1,teDerColIndices{2}(hourDiff)};
            teNow = measurements{i,teDerColIndices{2}(hourDiff)};
            if (~isempty(pPrev) && ~isempty(pNow))        
                measurements{i,pDerColIndices{3}(hourDiff)} = (pNow-pPrev) / double(hourDiff);
            end
            if (~isempty(tePrev) && ~isempty(teNow))                
                measurements{i,teDerColIndices{3}(hourDiff)} = (teNow-tePrev) / double(hourDiff);
            end        
        end
    end

    %% 4th order derivatives 
    for i=5:size(measurements,1)-4    
        i
        for j=1:size(hourDiffsMatter,2)
            hourDiff = hourDiffsMatter(j);
            pPrev = measurements{i-1,pDerColIndices{3}(hourDiff)};
            pNow = measurements{i,pDerColIndices{3}(hourDiff)};
            tePrev = measurements{i-1,teDerColIndices{3}(hourDiff)};
            teNow = measurements{i,teDerColIndices{3}(hourDiff)};
            if (~isempty(pPrev) && ~isempty(pNow))        
                measurements{i,pDerColIndices{4}(hourDiff)} = (pNow-pPrev) / double(hourDiff);
            end
            if (~isempty(tePrev) && ~isempty(teNow))                
                measurements{i,teDerColIndices{4}(hourDiff)} = (teNow-tePrev) / double(hourDiff);
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
                    k = k+1
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

        mkdir(fullfile(ISDStatisticsDestDir,stationName));
        stationFile = fopen(fullfile(ISDStatisticsDestDir,stationName,strcat(hour{1},'.txt')),'w');
        for k=1:size(values_new,2)
            fprintf(stationFile,'%s',[char(values_new{k}),char(newline)]);
        end
        fclose(stationFile);   
    end        
end