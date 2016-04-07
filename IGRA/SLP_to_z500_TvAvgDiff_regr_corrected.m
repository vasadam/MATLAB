function [z500, z500_e] = SLP_to_z500_TvAvgDiff_regr_corrected(SLP,...              % surface-level pressure (mbar)
                                                               elev,...             % elevation of the station
                                                               Tv,...               % virtual temperature (°K)
                                                               TvAvgDiff,...        % difference between surface-level Tv and average Tv (°K)
                                                               E,...                % standard error of TvAvgDiff
                                                               regrParamName,...    % name of the independent variable 'x' of the linear regression (Error = ax+b)
                                                               regrParamValue,...   % value of the independent variable 'x'
                                                               a,...                % 'a' regression coefficient (ax+b)
                                                               b,...                % 'b' regression coefficient
                                                               sigmaA,...           % sigma of 'a'
                                                               sigmaB)              % sigma of 'b'     
%SLP_to_z500_TvAvgDiff_regr_corrected Converts atmospheric pressure to z500 geopotential height

R = 287.058;                        % Specific gas constant of air, J/(kg*K)
g = 9.80665;                        % Gravitational acceleration, m/s^2    
ISAPressure = 1013.25;              % pressure of International Standard Atmosphere at mean sea level (hPa)
ISAEquivalentTemperature = 288.15;  % temperature of International Standard Atmosphere at mean sea level (K)
        
z500 = elev + R/g * (Tv+TvAvgDiff) * log(SLP/500);

% Regression coefficients were calculated by subtracting ISA parameter
% values from the measured ones to make errors sigma(a) and sigma(b) less
% confusing; so perform the subtraction now, too.
if (strcmp(regrParamName,'P'))
    regrParamValue = regrParamValue - ISAPressure;
elseif (strcmp(regrParamName,'Te'))
    regrParamValue = regrParamValue - ISAEquivalentTemperature;
end

% delta_z500 is the likely error of our estimation that is correlated with
% the regression variable, so we can try to mitigate it using linear
% regression.
delta_z500 = a*regrParamValue + b;
z500 = z500 - delta_z500;

% Calculate standard error of z500
TvAvg_min = Tv+TvAvgDiff - E;                           
TvAvg_max = Tv+TvAvgDiff + E;                           
z500_min = elev + R/g * TvAvg_min * log(SLP/500);
z500_max = elev + R/g * TvAvg_max * log(SLP/500);
z500_e = (z500_max - z500_min) / 2;         
end
