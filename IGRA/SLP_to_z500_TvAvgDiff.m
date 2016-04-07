function [z500, z500_e] = SLP_to_z500_TvAvgDiff(SLP,...       % surface-level pressure (mbar)
                                                elev,...      % elevation of the station
                                                Tv,...        % virtual temperature (°K)
                                                TvAvgDiff,... % difference between surface-level Tv and average Tv (°K)
                                                E)            % standard error of TvAvgDiff
%SLP_to_z500_TvAvgDiff Converts atmospheric pressure to z500 geopotential height

R = 287.058;    % Specific gas constant of air, J/(kg*K)
g = 9.80665;    % Gravitational acceleration, m/s^2    

z500 = elev + R/g * (Tv+TvAvgDiff) * log(SLP/500);

% Calculate standard error of z500
TvAvg_min = Tv+TvAvgDiff - E;                           
TvAvg_max = Tv+TvAvgDiff + E;                           
z500_min = elev + R/g * TvAvg_min * log(SLP/500);
z500_max = elev + R/g * TvAvg_max * log(SLP/500);
z500_e = (z500_max - z500_min) / 2;         
end
