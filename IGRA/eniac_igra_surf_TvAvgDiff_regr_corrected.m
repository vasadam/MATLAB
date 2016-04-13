clear    %    Clear the memory.
numStationsEurope = 125;
months_of_seasons = containers.Map;
months_of_seasons('winter') = [12 1 2];
months_of_seasons('spring') = [3 4 5];
months_of_seasons('summer') = [6 7 8];
months_of_seasons('autumn') = [9 10 11];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Basic parameters that may be varied
daylen  = 1;   %  Forecast length in days.
DtHours = 0.1;   %  Timestep in hours.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Define the time variables
seclen = daylen*24*60*60;        %  Integration time (in seconds)
dt = DtHours*60*60;             %  Timestep in seconds
nt = seclen/dt;       %  Total number of time-steps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Define some constants (SI units)
g = 9.80665;              % Gravitational acceleration (m/s^2)
a = (4*10^7)/(2*pi);      % Radius of the Earth
Omega = 2*pi/(24*60*60);  % Angular velocity of Earth.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Details of polar stereographic grid
LLCornerLat = 39;                           % latitude of the lower left (x=0,y=0) corner of the grid
LLCornerLon = -2.2;                         % longitude of the lower left (x=0,y=0) corner of the grid
LLCornerR = (cosd(LLCornerLat)/(sind(LLCornerLat)+1))*2*a;
LLCornerX = LLCornerR*sind(LLCornerLon);    % x coordinate of the lower left corner on the polar stereographic map
LLCornerY = -LLCornerR*cosd(LLCornerLon);   % y coordinate of the lower left corner on the polar stereographic map
size_x  = 10;                               % Points in x direction
size_y  = 10;                               % Points in y direction
ds = 300000;                                % Grid step at North Pole (m)
Centerangle = 0;                            % Central angle of map

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute latitude and longitude
%%  on the polar stereographic grid.
r2d = 180/pi;  % Conversion factor: radians to degrees.
d2r = pi/180;  % Conversion factor: degrees to radians.

X = zeros(size_y,size_x);
Y = zeros(size_y,size_x);
LON = zeros(size_y,size_x);
LAT = zeros(size_y,size_x);
for ny=1:size_y
    for nx=1:size_x
        X(nx,ny) = LLCornerX + (nx-1)*ds;
        Y(nx,ny) = LLCornerY + (ny-1)*ds;
        r = sqrt(X(nx,ny)^2+Y(nx,ny)^2);

        theta = atan2(Y(nx,ny),X(nx,ny));
        lambda = theta + d2r*(90+Centerangle);
        if(lambda>pi)
            lambda = lambda - 2*pi;
        end
        LON(nx,ny) = lambda;                        %  Longitude (radians)
        LAT(nx,ny) = 2*((pi/4)-atan(r/(2*a)));      %  Latitude  (radians)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute Coriolis Parameter, Map Factor
%%  and parameter h = g*m^2/f used in the BVE
FCOR = zeros(size_y,size_x);
H = zeros(size_y,size_x);
for ny=1:size_y
    for nx=1:size_x
        map = 2 / (1+sin(LAT(nx,ny)));
        FCOR(nx,ny) = 2*Omega*sin(LAT(nx,ny));
        H(nx,ny) = g * map^2 / FCOR(nx,ny);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Read the station parameters (ID,lat,lon)
stationlistfile = fopen('C:\Users\EDMMVAS\Documents\NOAA\igra2-station-list-filtered_Europe.txt');
keySet = cell(1,numStationsEurope);
valueSet = cell(1,numStationsEurope);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Read the TvAvg diferrence (mean,error) values for each (lat,lon) coordinate pair
TvAvgDiffRootDir = 'C:\Users\EDMMVAS\Documents\NOAA\IGRA_TvAvgDiff_from_2000';
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
% Remove empty tail of keySet2 and valueSet2
keySet2 = keySet2(1,1:m);
valueSet2 = valueSet2(1,1:m);
statMap = containers.Map(keySet2, valueSet2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read the regression coefficients which represent the connection between
%% surface-level atmospheric parameters (P or Te) and the error of the calculated value of z500
regrRootDir = 'C:\Users\EDMMVAS\Documents\NOAA\ISD_stat_tend_corr_regr_PandTe_per_season';
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Read the initial and verification data (TODO)
rootdir = 'D:\NOAA\IGRA_parsed_full';
yeardirs = dir(rootdir);
% Pre-allocate forecast array based on the number of years
forecasts = repmat(struct('datetime',[], ...    % date and time of the forecast
                          'Z0',[], ...          % initial z500 height field
                          'Z24_f',[], ...       % forecast z500 height field
                          'Z24_a',[], ...       % analysis z500 height field
                          'MAE_BVE',[], ...     % Mean Absolute Error of the z500 forecast (BVE)
                          'MAE_pers',[], ...    % Mean Absolute Error of the z500 forecast (persistence)
                          'P0',[], ...          % initial SLP field
                          'Z0_SLP',[], ...      % calculated initial z500 height field (from SLP)
                          'dZ0',[], ...         % error of the SLP->z500 conversion for initial z500 height field
                          'Z0_SLP_diff',[], ... % mean of absolute differences between Z0 and Z0_SLP
                          'Z24_f_SLP',[], ...   % forecast z500 height field (from Z0_SLP)
                          'Z24_a_SLP',[], ...   % calculated analysis z500 height field (from SLP)
                          'dZ24',[], ...        % error of the SLP->z500 conversion for analysis z500 height field
                          'P24_f',[], ...       % forecast SLP field (from Z24->SLP conversion)
                          'dP24',[], ...        % error of the z500->SLP conversion
                          'P24_a',[], ...       % analysis SLP field
                          'MAE_BVE_SLP',[], ... % Mean Absolute Error of the z500 forecast (BVE) based on SLP->z500
                          'MAE_pers_SLP',[] ... % Mean Absolute Error of the z500 forecast (persistence) based on SLP->z500
                     ),2*366*(size(yeardirs,1)-2), 1);
num = 1;
for i=1:size(yeardirs,1)
    if (~(yeardirs(i).isdir) || strcmp(yeardirs(i).name,'.') || strcmp(yeardirs(i).name,'..') ...  % Skip any files, '.' and '..'
        || (str2double(yeardirs(i).name) < 2008) )                   % Skip years before 2008
        continue;
    end
    daydirs = dir(fullfile(rootdir,yeardirs(i).name));
    
    for j=1:size(daydirs,1)        
        if (~daydirs(j).isdir || strcmp(daydirs(j).name,'.') || strcmp(daydirs(j).name,'..'))  % Skip any files, '.' and '..'
            continue;
        end
        
        fprintf('%s %s\n',yeardirs(i).name,daydirs(j).name);
        currentdaydir = daydirs(j);
        currentdayhourfiles = dir(fullfile(rootdir,yeardirs(i).name,currentdaydir.name));  
        
        if (j<size(daydirs,1))
            nextyeardirname = yeardirs(i).name;
            nextdaydir = daydirs(j+1);
            nextdayhourfiles = dir(fullfile(rootdir,yeardirs(i).name,nextdaydir.name));
        % If the current day is 31 December, look for the files from next
        % year's 1 January
        elseif (strcmp(currentdaydir.name,'12.31'))     
            nextyeardirname = num2str(str2double(yeardirs(i).name)+1);
            if (~exist(fullfile(rootdir,nextyeardirname,'01.01'),'dir'))
                break;
            end
            nextdaydir = struct('name','01.01','date','','bytes',[],'isdir',true,'datenum',[]);            
            nextdayhourfiles = dir(fullfile(rootdir,nextyeardirname,nextdaydir.name));     
        % Else skip this day
        else
            continue;
        end              
        
        for k=1:size(currentdayhourfiles)
            if (strcmp(currentdayhourfiles(k).name,'.') ...
                || strcmp(currentdayhourfiles(k).name,'..'))  % Skip '.' and '..'
                continue;
            end
            currentdayhourfile = currentdayhourfiles(k);
            % Find the same time file in the next day folder
            nextdayhourfile = struct;
            for l=1:size(nextdayhourfiles)
                if (strcmp(nextdayhourfiles(l).name,'.') ...
                    || strcmp(nextdayhourfiles(l).name,'..'))  % Skip '.' and '..'
                    continue;
                end                      
                if (strcmp(currentdayhourfile.name, nextdayhourfiles(l).name))
                    nextdayhourfile = nextdayhourfiles(l);
                    break;
                end
            end
            % Skip this hour if there is no matching file in the next day
            % folder
            if (isempty(fieldnames(nextdayhourfile)))
                continue;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Read the current and the next day files and perform grid
            % interpolation
            Data0 = load(fullfile(rootdir,yeardirs(i).name,currentdaydir.name,currentdayhourfile.name));
            LAT0 = Data0(:,1);
            LON0 = Data0(:,2);
            Z0 = Data0(:,3);
            ELEV0 = Data0(:,4);
            P0 = Data0(:,5);
            TV0 = Data0(:,6);
            TE0 = Data0(:,7);
            
            Data24 = load(fullfile(rootdir,nextyeardirname,nextdaydir.name,nextdayhourfile.name));
            LAT24 = Data24(:,1);
            LON24 = Data24(:,2);
            Z24 = Data24(:,3);
            ELEV24 = Data24(:,4);
            P24 = Data24(:,5);
            TV24 = Data24(:,6);            
            TE24 = Data24(:,7);           
            
            % Calculate (x,y) coordinates from (lat,lon)
            R0 = (cosd(LAT0)./(sind(LAT0)+1))*2*a;
            X0 = R0.*sind(LON0);
            Y0 = -R0.*cosd(LON0);
            
            R24 = (cosd(LAT24)./(sind(LAT24)+1))*2*a;
            X24 = R24.*sind(LON24);
            Y24 = -R24.*cosd(LON24);            
            
            % Perform grid interpolation on P,z500 and e(SLP->z500)
            F = scatteredInterpolant(X0,Y0,P0,'natural','linear');
            forecasts(num).P0 = F(X,Y); 
%             P0 = F(X,Y);            
            F = scatteredInterpolant(X0,Y0,Z0,'natural','linear');
            forecasts(num).Z0 = F(X,Y); 
%             Z0 = F(X,Y);       

            F = scatteredInterpolant(X24,Y24,P24,'natural','linear');
            forecasts(num).P24_a = F(X,Y); 
%             P24 = F(X,Y);  
            F = scatteredInterpolant(X24,Y24,Z24,'natural','linear');
            forecasts(num).Z24_a = F(X,Y); 
%             Z24 = F(X,Y);     

            % If interpolation fails due to insufficient number of data
            % points, skip to the next time file.
            if (isempty(forecasts(num).P0) || ...
                isempty(forecasts(num).Z0) || ...
                isempty(forecasts(num).P24_a) || ...
                isempty(forecasts(num).Z24_a))
                continue;
            end            
            
            % SLP->z500 conversion and interpolation for intitial and
            % analysis z500 height fields
            Z0_SLP = zeros(numel(P0),1);
            dZ0 = zeros(numel(P0),1);
            Z24_a_SLP = zeros(numel(P24),1);
            dZ24 = zeros(numel(P24),1);            
            for m=1:numel(P0)
                key = [num2str(LAT0(m),'%.4f') ...
                       num2str(LON0(m),'%.4f') ...
                       currentdaydir.name(1:2) ...
                       currentdayhourfile.name(1:2)];
                if (isKey(statMap,key))
                    statStruct = statMap(key);

                    % If there is regression data available, use it, otherwise
                    % use the simple method of SLP->z500 calculation
                    if (isfield(statStruct,'paramName'))
                        if (strcmp(statStruct.paramName,'P'))
                            paramValue = P0(m);
                        elseif (strcmp(statStruct.paramName,'Te'))
                            paramValue = TE0(m);
                        end       
                        [Z0_SLP(m),dZ0(m)] = SLP_to_z500_TvAvgDiff_regr_corrected(P0(m),...
                                                                                  ELEV0(m),...
                                                                                  TV0(m),...
                                                                                  statStruct.TvAvgDiff,...
                                                                                  statStruct.E,...
                                                                                  statStruct.paramName,...
                                                                                  paramValue,...
                                                                                  statStruct.a,...
                                                                                  statStruct.b,...
                                                                                  statStruct.sigmaA,...
                                                                                  statStruct.sigmaB);
                    else
                        [Z0_SLP(m),dZ0(m)] = SLP_to_z500_TvAvgDiff(P0(m),...
                                                                   ELEV0(m),...
                                                                   TV0(m),...
                                                                   statStruct.TvAvgDiff,...
                                                                   statStruct.E);
                    end
                end
            end
            for m=1:numel(P24)
                key = [num2str(LAT24(m),'%.4f') ...
                       num2str(LON24(m),'%.4f') ...
                       nextdaydir.name(1:2) ...
                       nextdayhourfile.name(1:2)];
                if (isKey(statMap,key))
                    statStruct = statMap(key);

                    % If there is regression data available, use it, otherwise
                    % use the simple method of SLP->z500 calculation
                    if (isfield(statStruct,'paramName'))
                        if (strcmp(statStruct.paramName,'P'))
                            paramValue = P24(m);
                        elseif (strcmp(statStruct.paramName,'Te'))
                            paramValue = TE24(m);
                        end       
                        [Z24_a_SLP(m),dZ24(m)] = SLP_to_z500_TvAvgDiff_regr_corrected(P24(m),...
                                                                                      ELEV24(m),...
                                                                                      TV24(m),...
                                                                                      statStruct.TvAvgDiff,...
                                                                                      statStruct.E,...
                                                                                      statStruct.paramName,...
                                                                                      paramValue,...
                                                                                      statStruct.a,...
                                                                                      statStruct.b,...
                                                                                      statStruct.sigmaA,...
                                                                                      statStruct.sigmaB);
                    else
                        [Z24_a_SLP(m),dZ24(m)] = SLP_to_z500_TvAvgDiff(P24(m),...
                                                                       ELEV24(m),...
                                                                       TV24(m),...
                                                                       statStruct.TvAvgDiff,...
                                                                       statStruct.E);
                    end
                end
            end
            F = scatteredInterpolant(X0,Y0,Z0_SLP,'natural','linear');
            forecasts(num).Z0_SLP = F(X,Y); 
            forecasts(num).Z0_SLP_diff = mean(abs(forecasts(num).Z0(:)-forecasts(num).Z0_SLP(:)));
%             forecasts(num).dZ0 = error_natural_neighbor(X0,Y0,dZ0,X,Y);   
            F = scatteredInterpolant(X24,Y24,Z24_a_SLP,'natural','linear');
            forecasts(num).Z24_a_SLP = F(X,Y); 
%             forecasts(num).dZ24 = error_natural_neighbor(X24,Y24,dZ24,X,Y);            
            
            % Run the forecast model for measured z500 and SLP->z500
            [Z24_f,MAE_BVE,MAE_pers] = eniac_core_2(forecasts(num).Z0, ...
                                                    forecasts(num).Z24_a, ...
                                                    ds, dt, nt, FCOR, H);
            hour = strrep(currentdayhourfile.name,'.',':');
            hour = hour(1:5);
            forecasts(num).datetime = strcat(yeardirs(i).name,'.',currentdaydir.name,'.',{' '},hour);
            forecasts(num).Z24_f = Z24_f;
            forecasts(num).MAE_BVE = MAE_BVE;
            forecasts(num).MAE_pers = MAE_pers;
            
            [Z24_f_SLP,MAE_BVE_SLP,MAE_pers_SLP] = eniac_core_2(forecasts(num).Z0_SLP, ...
                                                                forecasts(num).Z24_a_SLP, ...
                                                                ds, dt, nt, FCOR, H);
            hour = strrep(currentdayhourfile.name,'.',':');
            hour = hour(1:5);
            forecasts(num).datetime = strcat(yeardirs(i).name,'.',currentdaydir.name,'.',{' '},hour);
            forecasts(num).Z24_f_SLP = Z24_f_SLP;
            forecasts(num).MAE_BVE_SLP = MAE_BVE_SLP;
            forecasts(num).MAE_pers_SLP = MAE_pers_SLP;            
            num = num+1;
        end
    end
end

% Remove empty array elements
forecasts(num:size(forecasts)) = [];

