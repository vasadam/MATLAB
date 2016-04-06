function [z500, z500_e] = SLP_to_z500_grad(P,...     % pressure (mbar)
                                           elev,...  % elevation of the station
                                           Tv,...    % virtual temperature (°K)
                                           dT,...    % vertical temperature gradient (°K/m)
                                           dT_e)     % standard error of dT
%P_to_z500 Converts atmospheric pressure to z500 geopotential height

R = 287.058;    % Specific gas constant of air, J/(kg*K)
g = 9.80665;    % Gravitational acceleration, m/s^2        
% A = 17.625;
% B = 243.04;

C = (R/g)*log(P/500);
z500 = elev + C*(dT*elev-2*Tv)/(C*dT-2);

% Calculate standard error of z500
dT_min = dT - dT_e;                           
dT_max = dT + dT_e;                           
z500_min = elev + C*(dT_min*elev-2*Tv)/(C*dT_min-2);
z500_max = elev + C*(dT_max*elev-2*Tv)/(C*dT_max-2);
z500_e = (z500_max - z500_min) / 2;         
end

