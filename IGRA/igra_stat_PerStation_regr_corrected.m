% clear    %    Clear the memory.
% nl = java.lang.System.getProperty('line.separator').char;
USE_REGR_CORRECTION = false;

if (USE_REGR_CORRECTION)
    fprintf('!!!!!!!! REGR_CORRECTION: ON !!!!!!!!!\n');
else
    fprintf('!!!!!!!! REGR_CORRECTION: OFF !!!!!!!!!\n');
end

newline = java.lang.System.getProperty('line.separator');
numStationsEurope = 125;
months_of_seasons = containers.Map;
months_of_seasons('winter') = [12 1 2];
months_of_seasons('spring') = [3 4 5];
months_of_seasons('summer') = [6 7 8];
months_of_seasons('autumn') = [9 10 11];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Read the station parameters (ID,lat,lon)
stationListFile = fopen('D:\NOAA\igra2-station-list-filtered_Europe.txt','r');
keySet = cell(1,numStationsEurope);
valueSet = cell(1,numStationsEurope);
nextline = fgetl(stationListFile);
i=0;
while ischar(nextline)
   i = i+1;
   values = textscan(nextline,'%s\t%s\t%s\t%s\n');
   keySet{i} = values{1}{1};
   valueSet{i} = struct('lat',values{2}{1},'lon',values{3}{1});
   nextline = fgetl(stationListFile);
end
fclose(stationListFile);
stationMap = containers.Map(keySet,valueSet);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Read the TvAvg diferrence (mean,error) values for each (lat,lon) coordinate pair
TvAvgDiffRootDir = 'D:\NOAA\IGRA_TvAvgDiff_from_2000';
numAllStatFiles = size(rdir([TvAvgDiffRootDir, '\**\*stat.txt']),1);
keySet2 = cell(1,numAllStatFiles);
valueSet2 = cell(1,numAllStatFiles);
m=0;
stationdirs = dir(TvAvgDiffRootDir);
for i=1:size(stationdirs,1)
    if (strcmp(stationdirs(i).name,'.') || strcmp(stationdirs(i).name,'..')... % Skip '.' and '..'
        || ~any(~cellfun('isempty',strfind(keys(stationMap),stationdirs(i).name))))      % Skip stations not included in the station list file
        continue;
    end    
    
    monthdirs = dir(fullfile(TvAvgDiffRootDir,stationdirs(i).name));
    for j=1:size(monthdirs,1)
        if (strcmp(monthdirs(j).name,'.') || strcmp(monthdirs(j).name,'..'))  % Skip '.' and '..'
            continue;
        end        
        
        hourFiles = dir(fullfile(TvAvgDiffRootDir,stationdirs(i).name,monthdirs(j).name));
        for k=1:size(hourFiles)
            if (isempty(strfind(hourFiles(k).name,'stat')))    % Skip files whose name don't contain 'stat'
                continue;
            end
            
            statfile = fopen(fullfile(TvAvgDiffRootDir,stationdirs(i).name,monthdirs(j).name,hourFiles(k).name),'r');
            values = textscan(fgetl(statfile),'%s\t%s\n');      
            fclose(statfile);            
            m = m+1;
            fileNameParts = strsplit(hourFiles(k).name,'_');
            keySet2{m} = [stationMap(stationdirs(i).name).lat ...
                          stationMap(stationdirs(i).name).lon ...
                          monthdirs(j).name ...
                          fileNameParts{1}];    
            valueSet2{m} = struct('TvAvgDiff',str2double(values{1}),'E',str2double(values{2}));
        end        
    end
end
% Remove empty tail of keySet2 and valueSet2
keySet2 = keySet2(1,1:m);
valueSet2 = valueSet2(1,1:m);
statMap = containers.Map(keySet2, valueSet2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read the regression coefficients which represent the connection between
%% surface-level atmospheric parameters (P or Te) and the error of the calculated value of z500
regrRootDir = 'D:\NOAA\ISD_stat_tend_corr_regr_PandTe_per_season';
hourFiles = dir(fullfile(regrRootDir));
for i=1:size(hourFiles)
    if (strcmp(hourFiles(i).name,'.') || strcmp(hourFiles(i).name,'..'))  % Skip '.' and '..'
        continue;
    end
    
    statHourFile = fopen(fullfile(regrRootDir,hourFiles(i).name),'r');
    nextline = fgetl(statHourFile);
    j=0;
    while ischar(nextline)
        j = j+1;
        values = textscan(nextline,'%s %s %s %s %s %s %s %s %s %s\n','Delimiter',',');
        stationName = char(values{1});
        season = char(values{2});
        paramName = char(values{3});    % name of the independent variable of the regression
%         day = values{4};
%         corr = values{5};
%         pValue = values{6};
        a = str2double(values{7});
        b = str2double(values{8});        
        sigmaA = str2double(values{9});
        sigmaB = str2double(values{10});
        
        % Append regression coefficients to statStructs
        months = months_of_seasons(season);
        fileNameParts = strsplit(hourFiles(i).name,'.');
        key1 = [stationMap(stationName).lat ...
                stationMap(stationName).lon ...
                sprintf('%02d',months(1)) ...
                fileNameParts{1}];
        key2 = [stationMap(stationName).lat ...
                stationMap(stationName).lon ...
                sprintf('%02d',months(2)) ...
                fileNameParts{1}];
        key3 = [stationMap(stationName).lat ...
                stationMap(stationName).lon ...
                sprintf('%02d',months(3)) ...
                fileNameParts{1}];            

        if isKey(statMap,key1)
            statStruct1 = statMap(key1);
            statStruct1 = struct('TvAvgDiff',statStruct1.TvAvgDiff,...
                                 'E',statStruct1.E,...
                                 'paramName',paramName,...
                                 'a',a,...
                                 'b',b,...
                                 'sigmaA',sigmaA,...
                                 'sigmaB',sigmaB);
            statMap(key1) = statStruct1;
        end
        if isKey(statMap,key2)
            statStruct2 = statMap(key2);
            statStruct2 = struct('TvAvgDiff',statStruct2.TvAvgDiff,...
                                 'E',statStruct2.E,...
                                 'paramName',paramName,...
                                 'a',a,...
                                 'b',b,...
                                 'sigmaA',sigmaA,...
                                 'sigmaB',sigmaB);
            statMap(key2) = statStruct2;
        end
        if isKey(statMap,key3)
            statStruct3 = statMap(key3);
            statStruct3 = struct('TvAvgDiff',statStruct3.TvAvgDiff,...
                                 'E',statStruct3.E,...
                                 'paramName',paramName,...
                                 'a',a,...
                                 'b',b,...
                                 'sigmaA',sigmaA,...
                                 'sigmaB',sigmaB);
            statMap(key3) = statStruct3;
        end        
        
        nextline = fgetl(statHourFile);
    end
    fclose(statHourFile);
end

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read statistics and insert calculated z500 values
keySet3 = cell(1,2*numStationsEurope+2);    % per hour (per station + overall)
valueSet3 = cell(1,2*numStationsEurope+2);
m = 0;
numRecordsAll0 = 0;
numRecordsAll12 = 0;
sumAll0 = 0;
sumAll12 = 0;
PerStationStatisticsRootDir = 'D:\NOAA\IGRA_PerStationStatistics';
stationdirs = dir(PerStationStatisticsRootDir);
for i=1:size(stationdirs,1)
    if (strcmp(stationdirs(i).name,'.') || strcmp(stationdirs(i).name,'..')...  % Skip '.' and '..'
        || ~any(~cellfun('isempty',strfind(keys(stationMap),stationdirs(i).name))))       % Skip stations not included in the station list file
        continue;
    end        
    fprintf('%s\n',stationdirs(i).name);
    hourFiles = dir(fullfile(PerStationStatisticsRootDir,stationdirs(i).name));
    for j=1:size(hourFiles)
        if (strcmp(hourFiles(j).name,'.') || strcmp(hourFiles(j).name,'..'))  % Skip '.' and '..'
            continue;
        end        
               
        values_new = cell(1);        
        stationPos = stationMap(stationdirs(i).name);
        lat = stationPos.lat;
        lon = stationPos.lon;
        hourFile = fullfile(PerStationStatisticsRootDir,stationdirs(i).name,hourFiles(j).name);
        statfile = fopen(hourFile,'r');
        line = fgetl(statfile);
        k = 0;
        sum = 0;
        while(ischar(line) && size(line,2)>0)
            k = k+1;
            values = textscan(line,'%s %s %s %s %s %s %s %s %s %s %s %s\n');
        %     year = str2double(values{1});
            month = str2double(values{2});
        %     day = str2double(values{3});
            zLL = str2double(values{4});
            pLL = str2double(values{5});
        %     tLL = str2double(values{6});
        %     rhLL = str2double(values{7});
            tvLL = str2double(values{8});
            teLL = str2double(values{9});
            z500 = str2double(values{10});   
%             old_500_calc = str2double(values{11});
%             old_z500_calc_diff = str2double(values{12});

            fileNameParts = strsplit(hourFiles(j).name,'.');     
            key = [lat ...
                   lon ...
                   sprintf('%02d',month) ...
                   fileNameParts{1}];  
            if (isKey(statMap,key))
                statStruct = statMap(key);
                
                % If there is regression data available, use it, otherwise
                % use the simple method of SLP->z500 calculation
                if (USE_REGR_CORRECTION...
                    && isfield(statStruct,'paramName'))
                
                    if (strcmp(statStruct.paramName,'P'))
                        paramValue = pLL;
                    elseif (strcmp(statStruct.paramName,'Te'))
                        paramValue = teLL;
                    end
                    z500_calc = SLP_to_z500_TvAvgDiff_regr_corrected(pLL,...
                                                                     zLL,...
                                                                     tvLL,...
                                                                     statStruct.TvAvgDiff,...
                                                                     statStruct.E,...
                                                                     statStruct.paramName,...
                                                                     paramValue,...
                                                                     statStruct.a,...
                                                                     statStruct.b,...
                                                                     statStruct.sigmaA,...
                                                                     statStruct.sigmaB);                    
                else
                    z500_calc = SLP_to_z500_TvAvgDiff(pLL,...
                                                      zLL,...
                                                      tvLL,...
                                                      statStruct.TvAvgDiff,...
                                                      statStruct.E);
                end
 
                values_new{k} = cellstr([char(values{1}),' ',...
                                         char(values{2}),' ',...
                                         char(values{3}),' ',...
                                         char(values{4}),' ',...
                                         char(values{5}),' ',...
                                         char(values{6}),' ',...
                                         char(values{7}),' ',...
                                         char(values{8}),' ',...
                                         char(values{9}),' ',...
                                         char(values{10}),' ',...
                                         num2str(z500_calc,'%.0f'),' ',...
                                         num2str(z500_calc-z500,'%.0f')]);
                                     
                sum = sum + abs(z500_calc-z500);
                if (strcmp(fileNameParts{1},'00'))
                    sumAll0 = sumAll0 + abs(z500_calc-z500);
                elseif (strcmp(fileNameParts{1},'12'))
                    sumAll12 = sumAll12 + abs(z500_calc-z500);
                end
            else
                k = k-1; % Overwrite this line with the next one in the next iteration.
            end
            line = fgetl(statfile);
        end
        fclose(statfile);

        if (~exist(strrep(fullfile(PerStationStatisticsRootDir,stationdirs(i).name),'IGRA_PerStationStatistics','estimation_results'),'dir'))
            mkdir(strrep(fullfile(PerStationStatisticsRootDir,stationdirs(i).name),'IGRA_PerStationStatistics','estimation_results'));
        end
        statfile = fopen(strrep(hourFile,'IGRA_PerStationStatistics','estimation_results'),'w');
        for k=1:size(values_new,2)
            fprintf(statfile,'%s',[char(values_new{k}),char(newline)]);
        end
        fclose(statfile);
        
        if (strcmp(fileNameParts{1},'00'))
            numRecordsAll0 = numRecordsAll0 + k;
        elseif (strcmp(fileNameParts{1},'12'))    
            numRecordsAll12 = numRecordsAll12 + k;
        end
        
        m = m+1;
        keySet3{m} = [stationdirs(i).name,fileNameParts{1}];
        valueSet3{m} = num2str(sum / k, '%.2f');        
    end        
end

m = m+1;
keySet3{m} = ['overall','00'];
valueSet3{m} = num2str(sumAll0 / numRecordsAll0, '%.2f');
m = m+1;
keySet3{m} = ['overall','12'];
valueSet3{m} = num2str(sumAll12 / numRecordsAll12, '%.2f');
keySet3 = keySet3(1:m);
valueSet3 = valueSet3(1:m);
