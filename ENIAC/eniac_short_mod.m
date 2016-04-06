
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
% DELs = [ 184E+3 30000 30000 30000];
Cang = [ -90 -70 -35 -85 ]; % Angle between Date-line and
%                 			  positive y axis (varies).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Set the parameters for the forecast (Ncase)

M  = NUMx(Ncase); % Points in x direction
N  = NUMy(Ncase); % Points in y direction
Xp = Xpol(Ncase); % Coords of North Pole
Yp = Ypol(Ncase); % Coords of North Pole
Ds = DELs(Ncase); % Grid step at NP (metres)
Centerangle = Cang(Ncase); % Central angle of map.

MN   = M*N;       %  Total number of grid points.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute latitude and longitude
%%  on the polar stereographic grid.

r2d = 180/pi;  % Conversion factor: radians to degrees.
d2r = pi/180;  % Conversion factor: degrees to radians.

for ny=1:N
for nx=1:M
  xx = (nx-Xp)*Ds;
  yy = (ny-Yp)*Ds;
  rr = sqrt(xx^2+yy^2);
  theta = atan2(yy,xx);
  lambda = theta + d2r*(90+Centerangle);
  if(lambda>pi) lambda = lambda - 2*pi; end
  phi = 2*((pi/4)-atan(rr/(2*a)));
  LON(nx,ny) = lambda;   %  Longitude (radians)
  LAT(nx,ny) = phi;      %  Latitude  (radians)
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute Coriolis Parameter and Map Factor
%%  and parameter h = g*m^2/f used in the BVE
%%  (For psi-equation, h = m^2 is used).

for ny=1:N
for nx=1:M
  lambda = LON(nx,ny);
  phi = LAT(nx,ny);
  map = 2 / (1+sin(phi));
  f = 2*Omega*sin(phi);
  MAP(nx,ny) = map;
  FCOR(nx,ny) = f;
  h(nx,ny) = grav * map^2 / f;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Read and plot the initial and
%%%   verification height data

%% File name tag.
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
z0 = z;
z24 = load(File2);

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
%%  (First time step only)
if ( n == 1 )
   % Second x-derivative of z
    d2dx2z(2:M-1,:) = (z(3:M,:)+z(1:M-2,:)-2*z(2:M-1,:))/(Ds^2);

    % Second y-derivative of z
    d2dy2z(:,2:N-1) = (z(:,3:N)+z(:,1:N-2)-2*z(:,2:N-1))/(Ds^2); 
    
    %%  Laplacian of height
    %%  xi = d2dx2z + d2dy2z;
    xi(2:M-1,2:N-1) = d2dx2z(2:M-1,2:N-1)+d2dy2z(2:M-1,2:N-1);

    
    %%  Extend xi to boundaries (to compute Jacobian). xi from BCs in later steps
    xi(1,:) = 2*xi(2  ,:)-xi(3  ,:);  %  West
    xi(M,:) = 2*xi(M-1,:)-xi(M-2,:);  %  East
    xi(:,1) = 2*xi(:,2  )-xi(:,3  );  %  South
    xi(:,N) = 2*xi(:,N-1)-xi(:,N-2);  %  North
   
    %%  Copy initial vorticity field
    xinm1 = xi; 
%     xi0 = xi;
end

%%  Absolute vorticity
eta = h.*xi + FCOR;

% x-derivative of z
ddxz(2:M-1,:) = (z(3:M,:)-z(1:M-2,:))/(2*Ds);

% y-derivative of z
ddyz(:,2:N-1) = (z(:,3:N)-z(:,1:N-2))/(2*Ds);

% x-derivative of eta
ddxeta(2:M-1,:) = (eta(3:M,:)-eta(1:M-2,:))/(2*Ds);

% y-derivative of eta
ddyeta(:,2:N-1) = (eta(:,3:N)-eta(:,1:N-2))/(2*Ds);

%%  Compute the Jacobian J(eta,z)
ddtxi = ddxeta .* ddyz - ddyeta .* ddxz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Solve the Poisson Equation Del^2(ddtz) = ddtxi
%%  with homogeneous boundary conditions:  z is 
%%  constant on the boundaries, so ddtz vanishes.
% for num_iterations=1:(100*(M-2)*(N-2))    
%     i = randi([2,M-1]);
%     j = randi([2,N-1]);
%     ddtz(i,j) = (ddtz(i+1,j) + ddtz(i-1,j) + ddtz(i,j+1) + ddtz(i,j-1) - ddtxi(i,j)*(Ds^2))/4;     
% end

for num_it=1:100
    for i=2:(M-1)
        for j=2:(N-1)
            ddtz(i,j) = -(1/4)*(ddtxi(i,j)*(Ds^2) - ddtz(i+1,j) - ddtz(i-1,j) - ddtz(i,j+1) - ddtz(i,j-1));
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
sigW = max(0,+sign(ddyz(1,:)));
ddtxi(1,:) = sigW.*(2*ddtxi(2  ,:)-ddtxi(3  ,:)); 

%  Eastern boundary
sigE = max(0,-sign(ddyz(M,:)));
ddtxi(M,:) = sigE.*(2*ddtxi(M-1,:)-ddtxi(M-2,:));

%  Southern boundary
sigS = max(0,-sign(ddxz(:,1)));
ddtxi(:,1) = sigS.*(2*ddtxi(:,2  )-ddtxi(:,3  ));

%  Northern boundary
sigN = max(0,+sign(ddxz(:,N)));
ddtxi(:,N) = sigN.*(2*ddtxi(:,N-1)-ddtxi(:,N-2));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Step forward one time-step (leapfrog scheme).
xinp1 = xinm1 + (2*dt) * ddtxi;
znp1  = znm1  + (2*dt) * ddtz;

%  Move the old values into arrays znm1 and xinm1
%  Move the new values into arrays z and xi
znm1 = z;  
xinm1 = xi;  
z  = znp1;  
xi = xinp1;  

%  Restore the timestep (after the first step)
dt = Dt; 

if ( n==1 )
    z1 = z;
    eta1 = eta;
    xi1 = xinm1;
    ddtxi1 = ddtxi;
    ddtz1 = ddtz;
end

if (n==6)
    z6 = z;
end

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