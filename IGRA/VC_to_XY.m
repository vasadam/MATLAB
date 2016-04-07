function [X_vect,Y_vect] = VC_to_XY(V,C_vect)
% Converts Voronoi V and C_vect to X and Y vectors 

X_vect = zeros(1,numel(C_vect));
Y_vect = zeros(1,numel(C_vect));
for i=1:numel(C_vect)
    X_vect(i) = V(C_vect(i),1);
    Y_vect(i) = V(C_vect(i),2);
end
X_vect = {X_vect};
Y_vect = {Y_vect};
end

