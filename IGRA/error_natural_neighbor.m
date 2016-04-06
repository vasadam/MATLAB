function [E_interpolated] = error_natural_neighbor(X, Y, ...       % initial x and y coordinates (station coordinates)
                                                   E, ...          % initial error values
                                                   X_grid, Y_grid) % the grid to perform interpolation on
% Reshape the input points array
% Points = horzcat(reshape(X,size(X,1)*size(X,2),1),...
%                  reshape(Y,size(Y,1)*size(Y,2),1));
Points = horzcat(X,Y);

% Create initial voronoi diagram. We use two different versions of the
% command to eliminate the Inf values.
[V0,C0] = voronoin(Points);                     % contains cell information but also Inf values
[VX,VY] = voronoi(Points(:,1),Points(:,2));     % contains neither cell information nor Inf values
original_V0_size = size(V0,1);
for i=1:size(C0,1)
    [bool,index] = ismember(1,C0{i});
    if (bool)
        success = 0;
        neighbor_index = mod(index-2,size(C0{i},2))+1;
        neighbor_x = V0(C0{i}(neighbor_index),1);
        neighbor_y = V0(C0{i}(neighbor_index),2);
        for j=1:size(VX,2)
            % Find the neighbor point in VX,VY
            if (almost_equals(VX(1,j),neighbor_x) ...
                && almost_equals(VY(1,j),neighbor_y))
                    other_endpoint_x = VX(2,j);
                    other_endpoint_y = VY(2,j);
            elseif (almost_equals(VX(2,j),neighbor_x) ...
                    && almost_equals(VY(2,j),neighbor_y))
                        other_endpoint_x = VX(1,j);
                        other_endpoint_y = VY(1,j);   
            else
                continue;
            end   
            % If there is a match, check whether V0 contains the other
            % endpoint.
            other_endpoint_found = 0;
            for k=1:size(V0,1)
                if (almost_equals(V0(k,1),other_endpoint_x) ...
                    && almost_equals(V0(k,2),other_endpoint_y))
                        other_endpoint_found = 1;
                        break;
                end
            end
            % If V0 contains the other endpoint, check whether it was
            % originally there or was inserted by us. In the first case
            % it means that it's not Inf so that is not the one we're
            % looking for.
            
            % If V0 doesn't contain the other endpoint, then insert
            % other_endpoint_x and other_endpoint_y into V0
            if (~other_endpoint_found)
                V0 = vertcat(V0,[other_endpoint_x,other_endpoint_y]);
                k = size(V0,1);
            % If it was originally there, it's not Inf so continue.
            elseif (k<=original_V0_size)
                continue;
            end
            % Replace the "1" value in C0{i} with the inserted point's
            % index.
            C0{i}(index) = k;
            success = 1;
            break;
        end    
        i_str = num2str(i);
        sizev0str = num2str(size(V0,1));
        assert(success == 1, strcat('Failed to replace Inf value.',i_str,sizev0str));
    end
end

% For each grid point, calculate the natural neighbor interpolated
% error value
E_interpolated = zeros(size(X_grid,1),size(X_grid,2));
for i=1:size(X_grid,1)
    for j=1:size(X_grid,2)
        % Inf values are not likely to occur here because the grid points
        % should be inside the domain and not on the border (or outside).
        % So voronoi(...) is not necessary here.
        [V,C] = voronoin(vertcat(Points,[X_grid(i,j) Y_grid(i,j)]));

        % Find the intersecting cells with the grid point's cells and
        % store the intersection area values and weights for each.
        Values = [];
        Weights = [];
        [X_gpc,Y_gpc] = VC_to_XY(V,C{end});
        if (~ispolycw(X_gpc,Y_gpc))            
            [X_gpc,Y_gpc] = poly2cw(X_gpc,Y_gpc);
        end
        m = 0;
        for k=1:numel(C0)           
            [X_c,Y_c] = VC_to_XY(V0,C0{k});
            if (~ispolycw(X_c,Y_c)) 
                [X_c,Y_c] = poly2cw(X_c,Y_c);
            end
            [X_i,Y_i] = polybool('intersection', X_gpc,Y_gpc, X_c,Y_c);
            if (isempty(X_i))
                continue;
            end
            m = m+1;
            Values(m) = E(m);
            Weights(m) = polyarea(X_i{1},Y_i{1});
        end            
        Weights = Weights./sum(Weights);
        % Calculate error propagation from error values and weights
        E_interpolated(i,j) = sqrt(sum((Weights.^2).*(Values.^2)));
    end
end    
end

