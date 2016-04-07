function [z500, z500_e] = SLP_to_z500_multilevel(slp,...    % surface-level pressure (mbar)
                                                 elev,...   % elevation of the station
                                                 Tv_slp,... % surface-level virtual temperature (°K)
                                                 dTvs)      % array of dTv structs: {p,dTv,E}
%                                                  ps,...     % vector of pressures of available geopotential heights (mbar)
%                                                  dTvs,...   % vector of vertical virtual temperature gradients (°K/m)
%                                                  dTv_es)    % vector of standard errors of dTv
                                                 
%P_to_z500 Converts atmospheric pressure to z500 geopotential height

R = 287.058;    % Specific gas constant of air, J/(kg*K)
g = 9.80665;    % Gravitational acceleration, m/s^2        
% A = 17.625;
% B = 243.04;

Tv_prev = Tv_slp;
p_prev = slp;
z_prev = elev;
% for i=1:size(ps,2)
%     C = (R/g)*log(p_prev/ps(i));
%     z = z_prev + C*(dTvs(i)*z_prev-2*Tv_prev)/(C*dTvs(i)-2);      
%     if (i<size(ps,2))
%         Tv_prev = Tv_prev + dTvs(i)*(z-z_prev);
%         p_prev = ps(i);
%         z_prev = z;
%     end
% end
size(dTvs,2)
for i=1:size(dTvs,2)
    C = (R/g)*log(p_prev/dTvs(i).p);
    z = z_prev + C*(dTvs(i).dTv*z_prev - 2*Tv_prev) / (C*dTvs(i).dTv - 2);      
    if (i<size(dTvs,2))
        Tv_prev = Tv_prev + dTvs(i).dTv*(z-z_prev);
        p_prev = dTvs(i).p;
        z_prev = z;
    end
end
z500 = z;

% Calculate standard error of z500
% dTv_min = dTv - dT_e;                           
% dTv_max = dT + dT_e;                           
% z500_min = elev + C*(dTv_min*elev-2*Tv)/(C*dTv_min-2);
% z500_max = elev + C*(dTv_max*elev-2*Tv)/(C*dTv_max-2);
% z500_e = (z500_max - z500_min) / 2;         
z500_e = 0;
end

