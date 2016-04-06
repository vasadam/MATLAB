% clear    %    Clear the memory.
% nl = java.lang.System.getProperty('line.separator').char;
newline = java.lang.System.getProperty('line.separator');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Read the station parameters (ID,lat,lon)
stationlistfile = fopen('D:\NOAA\igra2-station-list-filtered.txt');
keySet = cell(1,694);
valueSet = cell(1,694);
nextline = fgetl(stationlistfile);
i=0;
while ischar(nextline)
   i = i+1;
   values = textscan(nextline,'%s\t%s\t%s\t%s\n');
   keySet{i} = values{1}{1};
   valueSet{i} = struct('lat',values{2}{1},'lon',values{3}{1});
   nextline = fgetl(stationlistfile);
end
fclose(stationlistfile);
stationMap = containers.Map(keySet,valueSet);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Read the TvAvg diferrence (mean,error) values for each (lat,lon) coordinate pair
TvAvgDiffRootDir = 'D:\NOAA\IGRA_TvAvgDiff_from_2000';
numStatFiles = size(rdir([TvAvgDiffRootDir, '\**\*stat.txt']),1);
keySet2 = cell(1,numStatFiles);
valueSet2 = cell(1,numStatFiles);
m=0;
stationdirs = dir(TvAvgDiffRootDir);
for i=1:size(stationdirs,1)
    if (strcmp(stationdirs(i).name,'.') || strcmp(stationdirs(i).name,'..'))  % Skip '.' and '..'
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
            
            statfile = fopen(fullfile(TvAvgDiffRootDir,stationdirs(i).name,monthdirs(j).name,hourFiles(k).name));
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
TvAvgDiffMap = containers.Map(keySet2, valueSet2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read the correlation and regression coefficients which represent the connection between
%% surface-level atmospheric parameters and the error of the calculated value of z500
corrRootDir = 'D:\NOAA\ISD_stat_tend_corr_per_season';
numCorrCoeffs = 4*size(rdir([corrRootDir, '\**\*.txt']),1); % take only the strongest correlation coeeficient for each station and season
keySet3 = cell(1,numCorrCoeffs);
valueSet3 = cell(1,numCorrCoeffs);
m=0;
stationdirs = dir(corrRootDir);
for i=1:size(stationdirs,1)
    if (strcmp(stationdirs(i).name,'.') || strcmp(stationdirs(i).name,'..'))  % Skip '.' and '..'
        continue;
    end       
        
    hourFiles = dir(fullfile(corrRootDir,stationdirs(i).name));
    for j=1:size(hourFiles)
        if (strcmp(hourFiles(j).name,'.') || strcmp(hourFiles(j).name,'..'))  % Skip '.' and '..'
            continue;
        end

        % Scan the file and find the strongest correlation for each season
        strongestCorrWinter = 0;
        strongestCorrSpring = 0;
        strongestCorrSummer = 0;
        strongestCorrAutumn = 0;
        corrFile = fopen(fullfile(corrRootDir,stationdirs(i).name,hourFiles(k).name));
        line = fgetl(corrFile);
        while(ischar(line))
            values = textscan(line,'%s,%s,%d,%f,%f\n');
        end
        values = textscan(fgetl(corrFile),'%s\t%s\n');      
        fclose(corrFile);
        
        m = m+1;
        fileNameParts = strsplit(hourFiles(j).name,'_');
        keySet3{m} = [stationMap(stationdirs(i).name).lat ...
                      stationMap(stationdirs(i).name).lon ...
                      monthdirs(j).name ...
                      fileNameParts{1}];    
        valueSet3{m} = struct('TvAvgDiff',str2double(values{1}),'E',str2double(values{2}));
    end        
end
CorrMap = containers.Map(keySet2, valueSet2);

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read statistics and insert calculated z500 values
PerStationStatisticsRootDir = 'D:\NOAA\IGRA_PerStationStatistics';
stationdirs = dir(PerStationStatisticsRootDir);
parfor i=1:size(stationdirs,1)
    if (strcmp(stationdirs(i).name,'.') || strcmp(stationdirs(i).name,'..'))  % Skip '.' and '..'
        continue;
    end        
    stationdirs(i).name
    hourFiles = dir(fullfile(PerStationStatisticsRootDir,stationdirs(i).name));
    for j=1:size(hourFiles)
        if (strcmp(hourFiles(j).name,'.') || strcmp(hourFiles(j).name,'..'))  % Skip '.' and '..'
            continue;
        end        

        values_new = cell(1);        
        stationPos = stationMap(stationdirs(i).name);
        lat = stationPos.lat;
        lon = stationPos.lon;
        hourFile = fullfile('D:\NOAA\IGRA_PerStationStatistics',stationdirs(i).name,hourFiles(j).name);
        statfile = fopen(hourFile,'r');
        line = fgetl(statfile);
        k = 0;
        while(ischar(line))
            k = k+1;
            values = textscan(line,'%s %s %s %s %s %s %s %s %s %s\n');
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

            fileNameParts = strsplit(hourFiles(j).name,'.');     
            key = [lat ...
                   lon ...
                   sprintf('%02d',month) ...
                   fileNameParts{1}];  
            if (isKey(TvAvgDiffMap,key))                
                TvAvgDiffStruct = TvAvgDiffMap(key);
                z500_calc = SLP_to_z500_TvAvgDiff(pLL, zLL, tvLL, TvAvgDiffStruct.TvAvgDiff, TvAvgDiffStruct.E);                
                values_new{k} = cellstr([line,' ',...
                                         num2str(z500_calc,'%.0f'),' ',...
                                         num2str(z500_calc-z500,'%.0f')]);                
            else
                k = k-1;    % Overwrite this line with the next one in the next iteration.
            end
            line = fgetl(statfile);
        end
        fclose(statfile);

        statfile = fopen(hourFile,'w');
        for k=1:size(values_new,2)
            fprintf(statfile,'%s',[char(values_new{k}),char(newline)]);
        end
        fclose(statfile);   
    end        
end
