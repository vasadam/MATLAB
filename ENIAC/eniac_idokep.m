
clear    %    Clear the memory.
clf      %    clear the display.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Basic parameters that may be varied

Ncase   = 1;   %  Forecast case number.
DAYLEN  = 1;   %  Forecast length in days.
DtHours = 1;   %  Timestep in hours.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Details of polar stereographic grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Set the parameters for the forecast (Ncase)

xmin = 420;           % x coordinate of the western border of the grid
xmax = 640;           % x coordinate of the eastern border of the grid
ymin = 180;           % y coordinate of the northern border of the grid
ymax = 400;           % y coordinate of the southern border of the grid
M  = 6;              % Points in x direction
N  = 6;              % Points in y direction
Ds = 30000;           % distance between grid points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Define the time variable

seclen = 900;        %  Integration time (in seconds)
Dt = 1;             %  Timestep in seconds
nt = seclen/Dt;       %  Total number of time-steps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Define some physical constants (SI units)

R = 287;            % Specific gas constant of air (J/(kg*K))
g = 9.80665;        % Gravitational acceleration (m/s^2)
coeff = R/g;        % Coefficient used for geopotential calculations
dT = 0.0065;        % Vertical temperature gradient (°C/m)
H500approx = 5750;  % Approximated 500 hPa geopotential height

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Set Coriolis Parameter
%%  and parameter h = g/f used in the BVE
f =  1.07e-04;      % Set Coriolis Parameter to a constant value for testing
for y = 1:N
    for x = 1:M
      FCOR(x,y) = f;
      h(x,y) = g / f;
    end
end

%%   Read the initial and verification data

P0 = load('idokep_P0.txt');
% P1 = load('idokep_P1.txt');
T0 = load('idokep_T0.txt');
% T1 = load('idokep_T1.txt');
RH0 = load('idokep_RH0.txt');
% RH1 = load('idokep_RH1.txt');

%%  Perform grid interpolation
[X,Y] = meshgrid(xmin:44:xmax,ymin:44:ymax);
F = scatteredInterpolant(P0(:,1),P0(:,2),P0(:,3),'natural','linear');
P0 = F(X,Y);
% F = scatteredInterpolant(P1(:,1),P1(:,2),P1(:,3),'natural','linear');
% P1 = F(X,Y);
F = scatteredInterpolant(T0(:,1),T0(:,2),T0(:,3),'natural','linear');
T0 = F(X,Y);
% F = scatteredInterpolant(T1(:,1),T1(:,2),T1(:,3),'natural','linear');
% T1 = F(X,Y);
F = scatteredInterpolant(RH0(:,1),RH0(:,2),RH0(:,3),'natural','linear');
RH0 = F(X,Y);
% F = scatteredInterpolant(RH1(:,1),RH1(:,2),RH1(:,3),'natural','linear');
% RH1 = F(X,Y);

%% Calculate the 500 hPa geopotential height
A = 17.625;         % coefficient A for dewpoint temperature calculation
B = 243.04;         % coefficient B for dewpoint temperature calculation (°C)
z = zeros(M,N);     % 500 hPa geopotential height
for x = 1 : size(X,1)
    for y = 1 : size(Y,1)
        Td = (B*(log(RH0(x,y)/100)+A*T0(x,y)/(B+T0(x,y)))) / (A-log(RH0(x,y)/100)-A*T0(x,y)/(B+T0(x,y)));
        Tv = (T0(x,y)+273.15) / (1-0.379*(6.11*(10^(7.5*Td/(237.7+Td)))/P0(x,y)));
        Tvavg = Tv-(H500approx*dT)/2;
        z(x,y) = coeff*Tvavg*log(P0(x,y)/500);
    end
end

%%  Read in the analyses on the grid
% z  = load(File1);  %  Initial height analysis
z0 = z;
% z24 = load(File2);

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

ztmp = zeros(60,M,N);

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
%%  (First time step only
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
    xi0 = xi;
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

ztmp(n,:,:) = z;

%  Restore the timestep (after the first step)
dt = Dt;
end     %%%%%    End of the time-stepping loop 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%     END MAIN LOOP    %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the (X,Y) grid (for plotting)
% [X, Y ] = meshgrid(1:M,1:N);
% X = X'; Y=Y';    % Transpose for standard plot.

% Plot the final height field
% zcontours = min(Z(:)):10:max(Z(:));
% contourf(X,Y,z,zcontours)
% title('FORECAST HEIGHT FIELD');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return