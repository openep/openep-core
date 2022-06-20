function [cv, cvX, interpCv] = doCvMapping_Triangulation( userdata)
% DOCVMAPPING_TRIANGULATION Calculates conduction velocities using
% triangulation of electrode activation times
%
% Usage:
%   [cv, cvX, interpCv] = doCvMapping_Triangulation( userdata, int )
% Where:
%   userdata - see importcarto_mem
%   int      - see openEpDataInterpolator.m
%   cv       - the calculated conduction velocity data, in m/s
%   cvX      - the Cartesian co-ordinates at which conduction velocity data
%              has been calculated. size(cvX) = [length(cv), 3].
%   interpCv - conduction velocity data interpolated across the surface of
%              the shell.
%              size(interpCv) = [length(userdata.surface.triRep.X), 1].
%
% DOCVMAPPING_TRIANGULATION Calculatess conduction velocities using
% triangulation
%
% Author: 
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% -----------------------------------s----------------------------
%
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% first perform global interpolation using triangulation to 
% calculate conduction velocities

% TODO

% accept only those conduction velocity values in proximity to electrodes?
% TODO. We should in some way limit cv and cvX only to values that are
% likely to be real; i.e. in close proximity; or at; mapping points.

% now do interpolation, using the interpolator specified, so we have a full
% dataset at each of the mesh nodes.
vtx = getVertices(userdata, 'used', false);

activationTime = getActivationTime(userdata);
activation_time_ms = activationTime * 1000;
[X, surfX] = getElectrogramX( userdata);
tri = getEgmTriangulation(userdata);

if (size(X,1) ~= size(tri.Points,1))
    error('dimensions does not match')
end

DT = delaunay(X(:,1:3));
DT1 = delaunayTriangulation(X);
TR = triangulation(DT1.ConnectivityList,X);
figure
hold on
trisurf(tri.ConnectivityList,tri.Points(:,1),tri.Points(:,2),tri.Points(:,3))
scatter3(X(:,1),X(:,2),X(:,3),20,'Color',[1 0 0])

neighbours_temp = [];

for elidx = 1:1:length(X)
    [n c v] = find(DT == elidx);
    for k = 1 :size(n,1)
        neighbours_temp = [neighbours_temp, DT(n(k),:)];
    end
    neighbours{elidx} = unique(neighbours_temp);
    neighbours_temp = [];
end
CV_total = [];
connected_electrodes = [];
CV_triangle = [];
for electrode_number = 1 : size(neighbours,2)
    center_electrode_index = ...
        find(neighbours{electrode_number} == electrode_number);
    center_electrode_number = neighbours{electrode_number}(center_electrode_index);
    
    for num_neighbours = 1 : size(neighbours{electrode_number},2)-2
        
        v1 = X(neighbours{electrode_number}(center_electrode_index),:);
        v2 = X(neighbours{electrode_number}(num_neighbours+1),:);
        v3 = X(neighbours{electrode_number}(num_neighbours+2),:);
        
        v31 = v3-v1;
        v21 = v2-v1;
        
        dist1 = sqrt(sum(v21.^2));
        dist2 = sqrt(sum(v31.^2));
        
        delta_T_21 = activation_time_ms(neighbours{electrode_number}(num_neighbours+1)) -  ...
            activation_time_ms(neighbours{electrode_number}(center_electrode_index));
        delta_T_31 = activation_time_ms(neighbours{electrode_number}(num_neighbours+2)) -  ...
            activation_time_ms(neighbours{electrode_number}(center_electrode_index));
        
        if dist1 <= 10 && dist1 >=1.5 && dist2 <= 10 && dist2 >= 1.5
            if delta_T_31 >= 2 && delta_T_21 >= 2 && (abs(dist1/delta_T_21)) >= 0.2 ...
                    && (abs(dist2/delta_T_31)) >= 0.2 && delta_T_31 < 30 && delta_T_21 < 30
                 %%%% Minimum activation differnce to be sure to
                 % not include electodes that are activated 
                 % simulatonously in the calculation
                   %
                theta_2 = atan2(norm(cross(v21,v31)),dot(v21,v31));
                if rad2deg(theta_2) >= 30
                    CV_triangle_temp = [center_electrode_number ...
                    neighbours{electrode_number}(num_neighbours+1) neighbours{electrode_number}(num_neighbours+2)];
                    CV_triangle = [CV_triangle; CV_triangle_temp];
                    quiver3(v1(1),v1(2),v1(3),v31(1),v31(2),v31(3),0)
                    quiver3(v1(1),v1(2),v1(3),v21(1),v21(2),v21(3),0)
                end
            end
        end
        
    end
end
CV =[];
CV_data = [];
for num_tri = 1 : size(CV_triangle,1)
    p = CV_triangle(num_tri,1);
    r = CV_triangle(num_tri,2);
    q = CV_triangle(num_tri,3);
    pq = sqrt(sum((X(p,:)- X(q,:)).^2));
    pr = sqrt(sum((X(p,:)- X(r,:)).^2));
    qr = sqrt(sum((X(q,:)- X(r,:)).^2));
    pq_v = X(q,:)- X(p,:);
    pr_v = X(r,:)- X(p,:);
    tpr = activation_time_ms(r)- activation_time_ms(p);
    tpq = activation_time_ms(q)- activation_time_ms(p);
    pr_vec = X(r,:)- X(p,:);
    pq_vec = X(q,:)- X(p,:);
    theta = acos((pr^2 + pq^2 - qr^2)/(2*pr*pq));
    alpha = atan((((tpq*pr)/(tpr*pq))-cos(theta))/sin(theta));
    CV_temp = (pq/tpq)*(cos(theta)*cos(alpha)+sin(theta)*sin(alpha));
    activation_time_p = [X(p,:), activation_time_ms(p)];
    activation_time_q = [X(q,:), activation_time_ms(q)];
    activation_time_r = [X(r,:), activation_time_ms(r)];
    
     if abs(CV_temp) > 0
         CV_data = [CV_data; [X(p,:),CV_temp]];
     end

end

Faces = userdata.surface.triRep.Triangulation;
Vertices = userdata.surface.triRep.X;
dist_list =[];
CV_list = [];

for num_vert = 1 : size(Vertices,1)
    for num_points = 1 : size(CV_data,1)
        dist = sqrt(sum((Vertices(num_vert,:)-CV_data(num_points,1:3)).^2));
        cv_temp_2 = CV_data(num_points,4);
        if dist < 10 && dist > 0
            dist_list = [dist_list dist];
            CV_list = [CV_list cv_temp_2];
        end
    end
    CV_vertice(num_vert) = mean(CV_list);
    dist_list = [];
    CV_list = [];
    
end
figure
hold on
ax1 = subplot(2,2,1);
patch(ax1,'Faces',Faces,'Vertices',Vertices,'FaceVertexCData',CV_vertice',...
       'FaceColor','interp','EdgeColor','none');

  
caxis([-0.2,2.0]);
colormap([[0.5,0.5,0.5];jet(10)]);
c2=colorbar;
ylim(c2,[0.0,2.0]);
ylabel(c2,'CV [m/s]');
axis off;
light;
lighting phong;
axis square
grid on;

ax2 = subplot(2,2,2);
histogram(ax2,CV_data(:,4),100)
xlim([0 2.2])
set(gca,'FontSize',20)
H=gca;
H.LineWidth=2;



ax3 = subplot(2,2,3);
hold on
patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',userdata.surface.act_bip(:,1),...
    'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
axis off;
light;
lighting phong;
axis square
grid on;
c3=colorbar;
ylabel(c3,'LAT (ms)');

end