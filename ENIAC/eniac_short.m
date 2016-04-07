
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Basic parameters that may be varied

    Ncase   = 1;   %  Forecast case number.
    DAYLEN  = 1;   %  Forecast length in days.
    DtHours = 1;   %  Timestep in hours.

%%  Indicator for psi-form of equation
    StreamFunction = 0;

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

t = (0:nt)*Dt;             %  Time variable
time = t/(60*60);          %  Time in hours (for plotting)
nspd = (24*60*60)/Dt;      %  Time steps per day

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Define some geophysical constants (SI units)

a = (4*10^7)/(2*pi);      %  Radius of the Earth
grav = 9.80665;           %  Gravitational acceleration
Omega = 2*pi/(24*60*60);  %  Angular velocity of Earth.
f0 = 2*Omega*sin(pi/4);   %  Mean Coriolis parameter.

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

  %%  Angles in degrees (for plotting)
  lamd = r2d*lambda; phid = r2d*phi;
  LONDEG(nx,ny) = lamd;  %  Longitude (degrees)
  LATDEG(nx,ny) = phid;  %  Latitude  (degrees)

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Record the latitude and longitude of the
%%  corners of the outer domain and of the
%%  inner domain (for checking purposes).
Corners = 0;
if(Corners == 1 )
  latSW = LATDEG(1,1); lonSW = LONDEG(1,1);
  latSE = LATDEG(M,1); lonSE = LONDEG(M,1);
  latNE = LATDEG(M,N); lonNE = LONDEG(M,N);
  latNW = LATDEG(1,N); lonNW = LONDEG(1,N);
          
  latSW = LATDEG(3  ,2  ); lonSW = LONDEG(3  ,2  );
  latSE = LATDEG(M-2,2  ); lonSE = LONDEG(M-2,2  );
  latNE = LATDEG(M-2,N-2); lonNE = LONDEG(M-2,N-2);
  latNW = LATDEG(3  ,N-2); lonNW = LONDEG(3  ,N-2);
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
  if (StreamFunction == 1) h(nx,ny) = map^2; end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Compute Sine Matrices for the Poisson solver

%%  Coefficients for x-transformation
for m1=1:M-2
for m2=1:M-2
  SM(m1,m2) = sin(pi*m1*m2/(M-1));
end
end

%%  Coefficients for y-transformation
for n1=1:N-2
for n2=1:N-2
  SN(n1,n2) = sin(pi*n1*n2/(N-1));
end
end

%%  Eigenvalues of Laplacian operator
for mm=1:M-2
for nn=1:N-2
  eigen = (sin(pi*mm/(2*(M-1))))^2 +(sin(pi*nn/(2*(N-1))))^2;
  EIGEN(mm,nn) = (-4/Ds^2) * eigen;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Read and plot the initial and
%%%   verification height data

%%%   Define the initial field sizes.
z0    = zeros(M,N);  %  Initial height
z24   = zeros(M,N);  %  Verifying Analysis
xi0   = zeros(M,N);  %  Initial Laplacian of height
eta0  = zeros(M,N);  %  Initial absolute vorticity

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
z0  = load(File1);  %  Initial height analysis
z24 = load(File2);  %  Verifying height analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Define working arrays to have correct size
ddxz    = zeros(M,N);     %  x derivative of z
ddyz    = zeros(M,N);     %  y derivative of z
gradzsq = zeros(M,N);     %  squared gradiant of z
d2dx2z  = zeros(M,N);     %  second x derivative of z
d2dy2z  = zeros(M,N);     %  second y derivative of z
xi      = zeros(M,N);     %  Laplacian of z
eta     = zeros(M,N);     %  absolute vorticity
ddxeta  = zeros(M,N);     %  x derivative of eta
ddyeta  = zeros(M,N);     %  y derivative of eta
Jacobi  = zeros(M,N);     %  Jacobian J(eta,z)
ddtz    = zeros(M,N);     %  Tendency of z
ddtxi   = zeros(M,N);     %  Tendency of xi
zdot    = zeros(M-2,N-2); %  Interior values of ddtz;
xidot   = zeros(M-2,N-2); %  Interior values of ddtxi;
ZDOT    = zeros(M-2,N-2); %  Fourier transform of zdot;
XIDOT   = zeros(M-2,N-2); %  Fourier transform of xidot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Install the initial conditions for the forecast

if(StreamFunction==1)
  %  z0    = grav * z0 ./ FCOR;  %  Initial streamfunction
  z0  = grav * z0  / f0;  %  Initial streamfunction
end

z = z0;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute the Laplacian of the height/psi
% Second x-derivative of z
d2dx2z(2:M-1,:) = (z(3:M,:)+z(1:M-2,:)-2*z(2:M-1,:))/(Ds^2);
% Second y-derivative of z
d2dy2z(:,2:N-1) = (z(:,3:N)+z(:,1:N-2)-2*z(:,2:N-1))/(Ds^2);

%%  Laplacian of height (or vorticity)
xi0(2:M-1,2:N-1) = d2dx2z(2:M-1,2:N-1)+d2dy2z(2:M-1,2:N-1);

%%  Extend xi0 to boundaries (for Jacobian).
xi0(1,:) = 2*xi0(2  ,:)-xi0(3  ,:);  %  West
xi0(M,:) = 2*xi0(M-1,:)-xi0(M-2,:);  %  East
xi0(:,1) = 2*xi0(:,2  )-xi0(:,3  );  %  South
xi0(:,N) = 2*xi0(:,N-1)-xi0(:,N-2);  %  North

%%  Absolute vorticity
eta0 = h.*xi0 + FCOR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%      MAIN LOOP     %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Integrate the BVE in time            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   Time-stepping is by leapfrog method
%%%   (first step is forward)
%%%
%%%   Define xi = Del^2(z). The BVE is
%%%
%%%       (d/dt)xi = J( h*xi+f, z)
%%%
%%%   We approximate the time derivative by
%%%     ( xi(n+1)-xi(n-1) ) / (2*Dt)    
%%%
%%%   The Jacobian term J(eta,z) is approximated
%%%   by centered space differences
%%%
%%%   Then the value of xi at the new time (n+1)*Dt is:
%%%      xi(n+1) =  xi(n-1) + 2*Dt * J(n)   
%%%
%%%   When we have ddtxi, we have to solve a Poisson
%%%   equation to get ddtz. Then both xi and z are stepped
%%%   forward to (n+1)*Dt and the cycle is repeated.
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Start of the time-stepping loop  %%%%%%%%%

for n=1:nt
      
%%  First time through the loop:
   if(n==1)
     dt = Dt/2;       %  First step is forward
     znm1  = z0;      %  Copy initial height field
     xinm1 = xi0;     %  Copy initial vorticity field
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute the derivatives, Laplacian and Jacobian.

% x-derivative of z
ddxz(2:M-1,:) = (z(3:M,:)-z(1:M-2,:))/(2*Ds);

% y-derivative of z
ddyz(:,2:N-1) = (z(:,3:N)-z(:,1:N-2))/(2*Ds);

% Square of the gradient of z
gradzsq = ddxz.^2+ddyz.^2;

% Second x-derivative of z
d2dx2z(2:M-1,:) = (z(3:M,:)+z(1:M-2,:)-2*z(2:M-1,:))/(Ds^2);

% Second y-derivative of z
d2dy2z(:,2:N-1) = (z(:,3:N)+z(:,1:N-2)-2*z(:,2:N-1))/(Ds^2);

%%  Laplacian of height
%%  xi = d2dx2z + d2dy2z;
xi(2:M-1,2:N-1) = d2dx2z(2:M-1,2:N-1)+d2dy2z(2:M-1,2:N-1);

%%  Extend xi to boundaries (to compute Jacobian).
%%  (First time step only; xi from BCs after that)
if ( n == 1 )
   xi(1,:) = 2*xi(2  ,:)-xi(3  ,:);  %  West
   xi(M,:) = 2*xi(M-1,:)-xi(M-2,:);  %  East
   xi(:,1) = 2*xi(:,2  )-xi(:,3  );  %  South
   xi(:,N) = 2*xi(:,N-1)-xi(:,N-2);  %  North
end

%%  Absolute vorticity
eta = h.*xi + FCOR;

% x-derivative of eta
ddxeta(2:M-1,:) = (eta(3:M,:)-eta(1:M-2,:))/(2*Ds);

% y-derivative of eta
ddyeta(:,2:N-1) = (eta(:,3:N)-eta(:,1:N-2))/(2*Ds);

%%  Compute the Jacobian J(eta,z)
Jacobi = ddxeta .* ddyz - ddyeta .* ddxz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculate the energy and enstrophy integrals
   E(n) = 0.5 * sum(sum(gradzsq));
   S(n) = 0.5 * sum(sum(xi.^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Solve the Poisson Equation Del^2(ddtz) = ddtxi
%%  with homogeneous boundary conditions:  z is 
%%  constant on the boundaries, so ddtz vanishes.

%% Note: Fourier transform of xidot denoted XIDOT
%%       Fourier transform of zdot  denoted ZDOT.
%%       Forward Fourier transform: 
%%          XIDOT = SM*xidot*SN;
%%       Inverse transform: 
%%          zdot = (4/((M-1)*(N-1)))*SM*ZDOT*SN;
%%

%  Tendency values in interior.
   xidot = Jacobi(2:M-1,2:N-1); 

%  Compute the transform of the solution
   XIDOT = SM*xidot*SN;

%  Convert transform of d(xi)/dt to transform of d(z)/dt
   ZDOT = XIDOT ./ EIGEN;
 
%  Compute inverse transform to get the height tendency.
   zdot = (4/((M-1)*(N-1))) * SM*ZDOT*SN;

   ddtz (2:M-1,2:N-1) = zdot;    %  Insert inner values 
   ddtxi(2:M-1,2:N-1) = xidot;   %  Insert inner values 

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Save the fields at quarterpoints of the integration
   if(n==1*nt/4) zq1 = z; end
   if(n==2*nt/4) zq2 = z; end
   if(n==3*nt/4) zq3 = z; end
   if(n==4*nt/4) zq4 = z; end

%  Restore the timestep (after the first step)
   dt = Dt; 

end     %%%%%    End of the time-stepping loop 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%     END MAIN LOOP    %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculate the energy and enstrophy
ddxz(2:M-1,:) = (z(3:M,:)-z(1:M-2,:))/(2*Ds);
ddyz(:,2:N-1) = (z(:,3:N)-z(:,1:N-2))/(2*Ds);
gradzsq = ddxz.^2+ddyz.^2;
E(nt+1) = 0.5 * sum(sum(gradzsq));
S(nt+1) = 0.5 * sum(sum(xi.^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert back from psi to z if necessary
if(StreamFunction==1)
  z0   = f0 * z0  / grav;  %  Initial streamfunction
  z    = f0 * z   / grav;  %  Initial streamfunction
  xi   = f0 * xi  / grav;  %  Initial Laplacian
  zq1  = f0 * zq1 / grav;  %  Solution at quarter-point
  zq2  = f0 * zq2 / grav;  %  Solution at quarter-point
  zq3  = f0 * zq3 / grav;  %  Solution at quarter-point
  zq4  = f0 * zq4 / grav;  %  Solution at quarter-point
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return