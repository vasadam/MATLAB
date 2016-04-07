% clear    %    Clear the memory.
% clf      %    clear the display.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Basic parameters that may be varied
Ncase   = 1;   %  Forecast case number.
DAYLEN  = 1;   %  Forecast length in days.
DtHours = 1;   %  Timestep in hours.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Details of polar stereographic grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Set the PS grid parameters for 4 cases.
NUMx = [ 19 19 19 19 ];   % points in x-direction
NUMy = [ 16 16 16 16 ];   % points in y-direction
Xpol = [ 10 10 10 10 ];   % x-coordinate of North Pole
Ypol = [ 14 12 14 14 ];   % y-coordinate of North Pole
DELs = [ 736E+3 736E+3 736E+3 736E+3];
Cang = [ -90 -70 -35 -85 ]; % Angle between Date-line and
                 			% positive y axis (varies).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Set the parameters for the forecast (Ncase)
M  = NUMx(Ncase); % Points in x direction
N  = NUMy(Ncase); % Points in y direction
Xp = Xpol(Ncase); % Coords of North Pole
Yp = Ypol(Ncase); % Coords of North Pole
Ds = DELs(Ncase); % Grid step at NP (metres)
Centerangle = Cang(Ncase); % Central angle of map.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Define the time variable
daylen = DAYLEN;           %  Integration time (in days)
seclen = daylen*24*60*60;  %  Integration time (in seconds)
Dt = DtHours*60*60;        %  Timestep in seconds
nt = seclen/Dt;            %  Total number of time-steps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Define some geophysical constants (SI units)
a = (4*10^7)/(2*pi);      %  Radius of the Earth
grav = 9.80665;           %  Gravitational acceleration
Omega = 2*pi/(24*60*60);  %  Angular velocity of Earth.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Read and plot the initial and
%%   verification height data
Case(1,1:6) = 'Case1-';
Case(2,1:6) = 'Case2-';
Case(3,1:6) = 'Case3-';
Case(4,1:6) = 'Case4-';

%% Initial analysis
YMDH1(1,1:10) = '1949010503';
YMDH1(2,1:10) = '1949013003';
YMDH1(3,1:10) = '1949013103';
YMDH1(4,1:10) = '1949021303';

%% Verifying analysis
YMDH2(1,1:10) = '1949010603';
YMDH2(2,1:10) = '1949013103';
YMDH2(3,1:10) = '1949020103';
YMDH2(4,1:10) = '1949021403';

%% Initial and verification analysis on PS grid
File1 = [Case(Ncase,:) YMDH1(Ncase,:) '.z00'];
File2 = [Case(Ncase,:) YMDH2(Ncase,:) '.z00'];

%%  Read in the analyses on the PS grid
z  = load(File1);  %  Initial height analysis
z24 = load(File2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute latitude and longitude
%%  on the polar stereographic grid.
%%  Compute fx, fy, hx, hy and h parameters
%%  used in the BVE.
r2d = 180/pi;  % Conversion factor: radians to degrees.
d2r = pi/180;  % Conversion factor: degrees to radians.

xxcoords = zeros(M,N);
yycoords = zeros(M,N);
fx = zeros(M,N);
fy = zeros(M,N);
hx = zeros(M,N);
hy = zeros(M,N);
h = zeros(M,N);
for ny=1:N
    for nx=1:M
      xx = (nx-Xp)*Ds;
      yy = (ny-Yp)*Ds;
      xxcoords(nx,ny) = xx;
      yycoords(nx,ny) = yy;
      rr = sqrt(xx^2+yy^2);
      theta = atan2(yy,xx);
      lambda = theta + d2r*(90+Centerangle);    % Longitude (radians)
      if (lambda > pi)
          lambda = lambda - 2*pi;
      end
      phi = 2*((pi/4)-atan(rr/(2*a)));          % Latitude  (radians)
      fx(nx,ny) = (-4*Omega*xx)/(1+xx^2+yy^2) * (1 + (1-xx^2-yy^2)/(1+xx^2+yy^2));
      fy(nx,ny) = (-4*Omega*yy)/(1+xx^2+yy^2) * (1 + (1-xx^2-yy^2)/(1+xx^2+yy^2));
      hx(nx,ny) = grav*xx/Omega * (1+xx^2+yy^2)/(1-xx^2-yy^2) * (2 + (1+xx^2+yy^2)/(1-xx^2-yy^2));
      hy(nx,ny) = grav*yy/Omega * (1+xx^2+yy^2)/(1-xx^2-yy^2) * (2 + (1+xx^2+yy^2)/(1-xx^2-yy^2));
      h(nx,ny) = grav/(2*Omega) * (1+xx^2+yy^2)^2/(1-xx^2-yy^2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute H matrices for each grid point
Hs = cell(M,N);
for i=2:(M-1)
    for j=2:(N-1)
        % Neighbors' coordinates
        xx_neighbors = [xxcoords(i-1,j), xxcoords(i-1,j+1), xxcoords(i,j+1), xxcoords(i+1,j+1),...
                        xxcoords(i+1,j), xxcoords(i+1,j-1), xxcoords(i,j-1), xxcoords(i-1,j-1)];
        yy_neighbors = [yycoords(i-1,j), yycoords(i-1,j+1), yycoords(i,j+1), yycoords(i+1,j+1),...
                        yycoords(i+1,j), yycoords(i+1,j-1), yycoords(i,j-1), yycoords(i-1,j-1)];
        H = zeros(8,8);
        for p=1:8
            H(p,:) = [xx_neighbors(p), yy_neighbors(p), xx_neighbors(p)^2, yy_neighbors(p)^2,...
                      xx_neighbors(p)^3, yy_neighbors(p)^3, xx_neighbors(p)^2 * yy_neighbors(p), xx_neighbors(p) * yy_neighbors(p)^2];
        end
        Hs{i,j} = H;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Define working arrays to have correct size
ddxz    = zeros(M,N);     %  x derivative of z
ddyz    = zeros(M,N);     %  y derivative of z
d2dx2z  = zeros(M,N);     %  second x derivative of z
d2dy2z  = zeros(M,N);     %  second y derivative of z
xi      = zeros(M,N);     %  Laplacian of z
eta     = zeros(M,N);     %  absolute vorticity
ddxeta  = zeros(M,N);     %  x derivative of eta
ddyeta  = zeros(M,N);     %  y derivative of eta
ddtz    = zeros(M,N);     %  Tendency of z
ddtxi   = zeros(M,N);     %  Tendency of xi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Install the initial conditions for the forecast
znm1 = z;        %  Copy initial height field
dt = Dt/2;       %  First step uses forward difference

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%      MAIN LOOP     %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Start of the time-stepping loop  %%%%%%%%%
for n=1:nt    
    %  Iterative algorithm to compute dz/dt    
    for num_it=1:100
        for i=2:(M-1)
            for j=2:(N-1)            
                z_neighbors = [z(i-1,j), z(i-1,j+1), z(i,j+1), z(i+1,j+1),...
                               z(i+1,j), z(i+1,j-1), z(i,j-1), z(i-1,j-1)]';
                a = Hs{i,j} \ z_neighbors;
                b = Hs{i,j} \ ones(8,1);
                c = a - b*z(i,j);
                ddtz(i,j) = 1/(-2*(b(3)+b(4))) * (2*(c(3)+c(4))*(hx(i,j)*c(2) - hy(i,j)*c(1))...
                                                  + (fx(i,j)*c(2) - fy(i,j)*c(1))...
                                                  + h(i,j)*((2*c(8) + 6*c(5))*c(2) - (2*c(7) + 6*c(6))*c(1)));                                          
            end
        end
    end  
    
    %  Step forward one time-step (leapfrog scheme).
    znp1  = znm1  + (2*dt) * ddtz;

    %  Move the old values into arrays znm1 and xinm1
    %  Move the new values into arrays z and xi
    znm1 = z;  
    z  = znp1;  

    %  Restore the timestep (after the first step)
    dt = Dt; 
end     %%%%%    End of the time-stepping loop 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%     END MAIN LOOP    %ff%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the (X,Y) grid (for plotting)
[X, Y ] = meshgrid(1:M,1:N);
X = X'; Y=Y';    % Transpose for standard plot.

% Plot the final height field
zcontours = 4500:50:6000;
contourf(X,Y,z,zcontours)
title('FORECAST HEIGHT FIELD');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return