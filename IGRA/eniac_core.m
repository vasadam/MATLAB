function [Z,MAE] = eniac_core(Z0, Z24, ...  % initial and analysis fields
                                      ds, ...        % Grid step at North Pole (m)
                                      dt, nt, ...    % timestep in seconds, total number of timesteps
                                      FCOR, H)       % Coriolis- and h parameters

% clear    %    Clear the memory.
% clf      %    clear the display.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Basic parameters that may be varied

% Ncase   = 1;   %  Forecast case number.
% daylen  = 1;   %  Forecast length in days.
% DtHours = 0.1; %  Timestep in hours.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Define some physical constants (SI units)

% R = 287;                  % Specific gas constant of air (J/(kg*K))
% g = 9.80665;              % Gravitational acceleration (m/s^2)
% coeff = R/g;              % Coefficient used for geopotential calculations
% dT = 0.0065;              % Vertical temperature gradient (°C/m)
% H500approx = 5750;        % Approximated 500 hPa geopotential height
% a = (4*10^7)/(2*pi);      % Radius of the Earth
% grav = 9.80665;           % Gravitational acceleration
% Omega = 2*pi/(24*60*60);  % Angular velocity of Earth.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Details of polar stereographic grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LLCornerLat = 39;  % latitude of the lower left (x=0,y=0) corner of the grid
% LLCornerLon = -2.2;  % longitude of the lower left (x=0,y=0) corner of the grid
% LLCornerR = (cosd(LLCornerLat)/(sind(LLCornerLat)+1))*2*a;
% LLCornerX = LLCornerR*sind(LLCornerLon);    % x coordinate of the lower left corner on the polar stereographic map
% LLCornerY = -LLCornerR*cosd(LLCornerLon);   % y coordinate of the lower left corner on the polar stereographic map
xsize  = size(Z0,2);    % Points in x direction
ysize  = size(Z0,1);    % Points in y direction
% ds = 300000;        % Grid step at North Pole (m)
% Centerangle = 0;    % Central angle of map

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute latitude and longitude
%%  on the polar stereographic grid.

% r2d = 180/pi;  % Conversion factor: radians to degrees.
% d2r = pi/180;  % Conversion factor: degrees to radians.
% 
% X = zeros(size_y,size_x);
% Y = zeros(size_y,size_x);
% LON = zeros(size_y,size_x);
% LAT = zeros(size_y,size_x);
% for ny=1:size_y
%     for nx=1:size_x
%         X(nx,ny) = LLCornerX + (nx-1)*ds;
%         Y(nx,ny) = LLCornerY + (ny-1)*ds;
%         r = sqrt(X(nx,ny)^2+Y(nx,ny)^2);
%         %   LON(nx,ny) = asin(xx(nx,ny)/rr(nx,ny));                                 %  Longitude (radians)
%         %   LAT(nx,ny) = acos((4*a/rr(nx,ny))*(rr(nx,ny)^2/(4*a^2+rr(nx,ny)^2)));   %  Latitude  (radians)
% 
%         theta = atan2(Y(nx,ny),X(nx,ny));
%         lambda = theta + d2r*(90+Centerangle);
%         if(lambda>pi) lambda = lambda - 2*pi; end
%         LON(nx,ny) = lambda;   %  Longitude (radians)
%         LAT(nx,ny) = 2*((pi/4)-atan(r/(2*a)));      %  Latitude  (radians)
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Define the time variable

% seclen = daylen*24*60*60;        %  Integration time (in seconds)
% dt = DtHours*60*60;             %  Timestep in seconds
% nt = seclen/dt;       %  Total number of time-steps.

%%   Read the initial and verification data

% Z0 = load('00.00.txt');
% LAT2 = Z0(:,1);
% LON2 = Z0(:,2);
% Z0 = Z0(:,3);
% R = (cosd(LAT2)./(sind(LAT2)+1))*2*a;
% X = R.*sind(LON2);
% Y = -R.*cosd(LON2);
% F = scatteredInterpolant(X,Y,Z0,'natural','linear');
% Z = F(X,Y);
% P0 = load('idokep_P0.txt');
% % P1 = load('idokep_P1.txt');
% T0 = load('idokep_T0.txt');
% % T1 = load('idokep_T1.txt');
% RH0 = load('idokep_RH0.txt');
% % RH1 = load('idokep_RH1.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute Coriolis Parameter and Map Factor
%%  and parameter h = g*m^2/f used in the BVE
%%  (For psi-equation, h = m^2 is used).

% FCOR = zeros(size_y,size_x);
% H = zeros(size_y,size_x);
% for ny=1:size_y
%     for nx=1:size_x
%         map = 2 / (1+sin(LAT(nx,ny)));
%         FCOR(nx,ny) = 2*Omega*sin(LAT(nx,ny));
%         H(nx,ny) = g * map^2 / FCOR(nx,ny);
%     end
% end

%%  Perform grid interpolation
% [X,Y] = meshgrid(xmin:44:xmax,ymin:44:ymax);
% F = scatteredInterpolant(P0(:,1),P0(:,2),P0(:,3),'natural','linear');
% P0 = F(X,Y);
% % F = scatteredInterpolant(P1(:,1),P1(:,2),P1(:,3),'natural','linear');
% % P1 = F(X,Y);
% F = scatteredInterpolant(T0(:,1),T0(:,2),T0(:,3),'natural','linear');
% T0 = F(X,Y);
% % F = scatteredInterpolant(T1(:,1),T1(:,2),T1(:,3),'natural','linear');
% % T1 = F(X,Y);
% F = scatteredInterpolant(RH0(:,1),RH0(:,2),RH0(:,3),'natural','linear');
% RH0 = F(X,Y);
% % F = scatteredInterpolant(RH1(:,1),RH1(:,2),RH1(:,3),'natural','linear');
% % RH1 = F(X,Y);

% %% Calculate the 500 hPa geopotential height
% A = 17.625;         % coefficient A for dewpoint temperature calculation
% B = 243.04;         % coefficient B for dewpoint temperature calculation (°C)
% z = zeros(M,N);     % 500 hPa geopotential height
% for x = 1 : size(X,1)
%     for y = 1 : size(Y,1)
%         Td = (B*(log(RH0(x,y)/100)+A*T0(x,y)/(B+T0(x,y)))) / (A-log(RH0(x,y)/100)-A*T0(x,y)/(B+T0(x,y)));
%         Tv = (T0(x,y)+273.15) / (1-0.379*(6.11*(10^(7.5*Td/(237.7+Td)))/P0(x,y)));
%         Tvavg = Tv-(H500approx*dT)/2;
%         z(x,y) = coeff*Tvavg*log(P0(x,y)/500);
%     end
% end

%%  Read in the analyses on the grid
% z  = load(File1);  %  Initial height analysis
% z0 = z;
% z24 = load(File2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Define working arrays to have correct size
DDXZ    = zeros(xsize,ysize);     %  x derivative of z
DDYZ    = zeros(xsize,ysize);     %  y derivative of z
D2DX2Z  = zeros(xsize,ysize);     %  second x derivative of z
D2DY2Z  = zeros(xsize,ysize);     %  second y derivative of z
XI      = zeros(xsize,ysize);     %  Laplacian of z
% ETA     = zeros(size_x,size_y);     %  absolute vorticity
DDXETA  = zeros(xsize,ysize);     %  x derivative of eta
DDYETA  = zeros(xsize,ysize);     %  y derivative of eta
DDTZ    = zeros(xsize,ysize);     %  Tendency of z
% DDTXI   = zeros(size_x,size_y);     %  Tendency of xi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Install the initial conditions for the forecast
Z = Z0;
ZNM1 = Z0;       %  Copy initial height field
dt = dt/2;       %  First step uses forward difference

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%      MAIN LOOP     %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Start of the time-stepping loop  %%%%%%%%%
for n=1:nt
%%  (First time step only
if ( n == 1 )
   % Second x-derivative of z
    D2DX2Z(2:xsize-1,:) = (Z(3:xsize,:)+Z(1:xsize-2,:)-2*Z(2:xsize-1,:))/(ds^2);

    % Second y-derivative of z
    D2DY2Z(:,2:ysize-1) = (Z(:,3:ysize)+Z(:,1:ysize-2)-2*Z(:,2:ysize-1))/(ds^2); 
    
    %%  Laplacian of height
    %%  xi = d2dx2z + d2dy2z;
    XI(2:xsize-1,2:ysize-1) = D2DX2Z(2:xsize-1,2:ysize-1)+D2DY2Z(2:xsize-1,2:ysize-1);
    
    %%  Extend xi to boundaries (to compute Jacobian). xi from BCs in later steps
    XI(1,:) = 2*XI(2,:)-XI(3 ,:);  %  West
    XI(xsize,:) = 2*XI(xsize-1,:)-XI(xsize-2,:);  %  East
    XI(:,1) = 2*XI(:,2)-XI(:,3);  %  South
    XI(:,ysize) = 2*XI(:,ysize-1)-XI(:,ysize-2);  %  North
   
    %%  Copy initial vorticity field
    XINM1 = XI; 
%     XI0 = XI;
end

%%  Absolute vorticity
ETA = H.*XI + FCOR;

% x-derivative of z
DDXZ(2:xsize-1,:) = (Z(3:xsize,:)-Z(1:xsize-2,:))/(2*ds);

% y-derivative of z
DDYZ(:,2:ysize-1) = (Z(:,3:ysize)-Z(:,1:ysize-2))/(2*ds);

% x-derivative of eta
DDXETA(2:xsize-1,:) = (ETA(3:xsize,:)-ETA(1:xsize-2,:))/(2*ds);

% y-derivative of eta
DDYETA(:,2:ysize-1) = (ETA(:,3:ysize)-ETA(:,1:ysize-2))/(2*ds);

%%  Compute the Jacobian J(eta,z)
DDTXI = DDXETA .* DDYZ - DDYETA .* DDXZ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Solve the Poisson Equation Del^2(ddtz) = ddtxi
%%  with homogeneous boundary conditions:  z is 
%%  constant on the boundaries, so ddtz vanishes.
for num_it=1:100
    for i=2:(xsize-1)
        for j=2:(ysize-1)
            DDTZ(i,j) = -(1/4)*(DDTXI(i,j)*(ds^2) - DDTZ(i+1,j) - DDTZ(i-1,j) - DDTZ(i,j+1) - DDTZ(i,j-1));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute ddtxi on the boundaries.
%%  If fluid is entering through the boundary,
%%  we set ddtxi to zero. If fluid is leaving 
%%  the region, we extrapolate ddtxi linearly
%%  from the interior. Charney, et al., Eqn (21).

%  Western boundary
sigW = max(0,+sign(DDYZ(1,:)));
DDTXI(1,:) = sigW.*(2*DDTXI(2  ,:)-DDTXI(3  ,:)); 

%  Eastern boundary
sigE = max(0,-sign(DDYZ(xsize,:)));
DDTXI(xsize,:) = sigE.*(2*DDTXI(xsize-1,:)-DDTXI(xsize-2,:));

%  Southern boundary
sigS = max(0,-sign(DDXZ(:,1)));
DDTXI(:,1) = sigS.*(2*DDTXI(:,2  )-DDTXI(:,3  ));

%  Northern boundary
sigN = max(0,+sign(DDXZ(:,ysize)));
DDTXI(:,ysize) = sigN.*(2*DDTXI(:,ysize-1)-DDTXI(:,ysize-2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Step forward one time-step (leapfrog scheme).
xinp1 = XINM1 + (2*dt) * DDTXI;
znp1  = ZNM1  + (2*dt) * DDTZ;

%  Move the old values into arrays znm1 and xinm1
%  Move the new values into arrays z and xi
ZNM1 = Z;  
XINM1 = XI;  
Z = znp1;  
XI = xinp1;  

%  Restore the timestep (after the first step)
if (n==1)
    dt = dt*2;
end

end     %%%%%    End of the time-stepping loop 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%     END MAIN LOOP    %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate the Forecast RMSE
ERROR = abs(Z24(2:ysize-1,2:xsize-1)-Z(2:ysize-1,2:xsize-1));
MAE = sum(ERROR(:)) / ((ysize-2)*(xsize-2));

% Define the (X,Y) grid (for plotting)
% [X, Y ] = meshgrid(1:M,1:N);
% X = X'; Y=Y';    % Transpose for standard plot.

% Plot the final height field
% zmin = floor(min(Z(:)));
% zmax = floor(max(Z(:)));
% zstep = 10;
% zcontours = zmin-mod(zmin,zstep) : zstep : zmax-mod(zmax,zstep)+zstep;
% contourf(X,Y,Z,zcontours)
% title('FORECAST HEIGHT FIELD');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%