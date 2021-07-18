function [locations, gradients] = trigrad(tri, sc)
% TRIGRAD calculates the gradient of a scalar field across a TriRep object
% Usage:
%   [locations gradients] = trigrad(tri, sc)
% Where:
%   locations contains the position at the centroid of each face
%   gradients contains the vector of grad(sc), corresponding to locations
%   tri is the trirep object containing the geometry
%   sc is the value of the scalar field at each point in tri.X
%   g is the gradient
% TRIGRAD calculates the gradient of sc across the surface defined by tri.
% The geometry can be a surface in 2d or 3d.

% Author: Nick Linton (2009)
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------
                        % %intiate random things
                        % if true
                        %     nPoints = 2000;
                        %     p = rand(nPoints,3)*2*pi;
                        %     for i = 1:nPoints
                        %         p(i,:) = [1 1 1] * rotationmatrix(p(i,1), p(i,2), p(i,3));
                        %     end
                        % end
                        % p(:,1) = p(:,1)-1;
                        % 
                        % tri = convhulln(p);
                        % tri = TriRep(tri,p);
                        % 
                        % SENSITIVITY = 0.1;
                        % 
                        % for i = 1:length(tri.X(:,1))
                        %     potential = tri.X(:,1) .* exp( -tri.X(:,1).^2 - tri.X(:,2).^2 );
                        % end
                        % [locations gradients] = trigrad(tri, potential);
                        % 
                        % close all
                        % figure
                        % hold on
                        % colormap(jet(64))
                        % c = (potential-min(potential)) ./ (max(potential)-min(potential)) * 63  +  1;
                        % h = trisurf(tri);
                        % set(h, 'CData',c ,'CDataMapping','direct', 'FaceColor', 'interp')
                        % for i = 1:size(locations,1)
                        %     start = locations(i,:);
                        %     loc = [start; start];
                        %     loc(2,:) = loc(2,:) + SENSITIVITY * gradients(i,:);
                        %     
                        %     plot3(start(:,1), start(:,2), start(:,3), 'Marker','*')
                        %     plot3(loc(:,1), loc(:,2), loc(:,3))
                        % end
                        % axis equal
                        % axis vis3d


    [~, locations] = tricentroid(tri);
    gradients = zeros(size(tri.Triangulation,1), 3);
    x = tri.X;
    nDim = size(x,2);
    if nDim == 2
        %then we only have 2 dimensions - so expand to 3, and retract later
        x = [x zeros(size(x,1),1)];
    end
    for i = 1:size(tri.Triangulation,1)
        gradients(i,:) = grad(x(tri.Triangulation(i,:),:), sc(tri.Triangulation(i,:)));
    end
    if nDim == 2
        gradients(:,3) = [];
    end
end

function g = grad(p, s)
    % first find two vectors that are orthogonal and on the surface of p1/2/3
    u = p(2,:)-p(1,:);
    u = u / norm(u);
    %u = reshape(u,1,length(u)); %make sure it's a row vector [it always is]

    v = p(3,:)-p(1,:);
    v = v - u.*dot(u,v);   %v is the component of v normal to u (u has already been normalised)
    v = v / norm(v);
    %v = reshape(v,1,length(v)); %make sure it's a row vector [it always is]

    % now express the change in s across the triangle 
    % ds = a1 * dp.u  +  a2 * dp.v
    a = [ dot(p(2,:)-p(1,:),u) , dot(p(2,:)-p(1,:),v) ;  dot(p(3,:)-p(1,:),u) , dot(p(3,:)-p(1,:),v) ] \ [ (s(2)-s(1)) ; (s(3)-s(1)) ];
    % so now, change in activation time is given by
    % s-s(1) = [u v]*a
    % so the vector [u v]*a is the gradient of s
    g = a(1)*u + a(2)*v;
end

            
            
            
            
            