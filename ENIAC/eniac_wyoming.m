clear    %    Clear the memory.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Basic parameters that may be varied
daylen  = 1;   %  Forecast length in days.
DtHours = 0.05;   %  Timestep in hours.

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
        %   LON(nx,ny) = asin(xx(nx,ny)/rr(nx,ny));                                 %  Longitude (radians)
        %   LAT(nx,ny) = acos((4*a/rr(nx,ny))*(rr(nx,ny)^2/(4*a^2+rr(nx,ny)^2)));   %  Latitude  (radians)

        theta = atan2(Y(nx,ny),X(nx,ny));
        lambda = theta + d2r*(90+Centerangle);
        if(lambda>pi) lambda = lambda - 2*pi; end
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
%%   Read the initial and verification data (TODO)
rootdir = 'D:\Wyoming_Z';
yeardirs = dir(rootdir);
% Pre-allocate forecast array based on the number of years
forecasts = repmat(struct('datetime',[],'Z',[],'MAE',[]),366*(size(yeardirs,1)-2), 1);
num = 1;
for i=1:size(yeardirs,1)
    if (~(yeardirs(i).isdir) || strcmp(yeardirs(i).name,'.'))  % Skip any files, '.' and '..'
        continue;
    end
    daydirs = dir(fullfile(rootdir,yeardirs(i).name));
    
    for j=1:size(daydirs,1)        
        if (~daydirs(j).isdir || strcmp(daydirs(j).name,'.') || strcmp(daydirs(j).name,'..'))  % Skip any files, '.' and '..'
            continue;
        end
        
        currentdaydir = daydirs(j)
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
            % Read the current and the next day files into Z0 and Z24
            Z0 = load(fullfile(rootdir,yeardirs(i).name,currentdaydir.name,currentdayhourfile.name));
            LAT0 = Z0(:,1);
            LON0 = Z0(:,2);
            Z0 = Z0(:,3);
            R0 = (cosd(LAT0)./(sind(LAT0)+1))*2*a;
            X0 = R0.*sind(LON0);
            Y0 = -R0.*cosd(LON0);
            F = scatteredInterpolant(X0,Y0,Z0,'natural','linear');
            Z0 = F(X,Y);           
            
            Z24 = load(fullfile(rootdir,nextyeardirname,nextdaydir.name,nextdayhourfile.name));
            LAT24 = Z24(:,1);
            LON24 = Z24(:,2);
            Z24 = Z24(:,3);
            R24 = (cosd(LAT24)./(sind(LAT24)+1))*2*a;
            X24 = R24.*sind(LON24);
            Y24 = -R24.*cosd(LON24);
            F = scatteredInterpolant(X24,Y24,Z24,'natural','linear');
            Z24 = F(X,Y);     
            
            % If interpolation fails due to insufficient number of data
            % points, skip to the next time file.
            if (isempty(Z0) || isempty(Z24))
                continue;
            end
            
            % Run the forecast model
            [Z,MAE] = eniac_wyoming_core(Z0, Z24, ds, dt, nt, FCOR, H);
            hour = strrep(currentdayhourfile.name,'.',':');
            hour = hour(1:5);
            forecasts(num).datetime = strcat(yeardirs(i).name,'.',currentdaydir.name,'.',{' '},hour);
            forecasts(num).Z = Z;
            forecasts(num).MAE = MAE;
            num = num+1;
        end
    end
end

% Remove empty array elements
forecasts(num:size(forecasts)) = [];
