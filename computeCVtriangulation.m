%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes CV magnitude and direction using triangulation method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v v_vec] = computeCVtriangulation(p,tp,q,tq,r,tr)
% inputs:
% coordinate vectors of 3 points defining a triangle: p, q, r i.e. p=[px,py,pz]
% local activation timings of each point: tp, tq, tr
% outputs CV as magntidue as well as normalised velocity vector

% Defines lengths of triange edges
x_pq = [q(1)-p(1);q(2)-p(2);q(3)-p(3)];
x_pr = [r(1)-p(1);r(2)-p(2);r(3)-p(3)];
x_qr = [r(1)-q(1);r(2)-q(2);r(3)-q(3)];

% Defines relative activation timings
t_pq = abs(tp - tq);
t_pr = abs(tp - tr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes velocity magnitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes angle theta (between x_pq and x_pr) using cosine rule
theta = (norm(x_pq)^2 + norm(x_pr)^2 - norm(x_qr)^2)/(2*norm(x_pq)*norm(x_pr));

% Computes angle alpha (between x_pq and velocity vector)
alpha = arctan( ( (t_pr*norm(x_pq))/(t_pq*norm(x_pr)) - cos(theta) ) / sin(theta) );

% Computes velocity
v = norm(x_pq)*cos(theta)/t_pq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes velocity as vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defines normal to plane (as unit vector)
n_plane_full = cross(x_pr,x_pq);
n_plane = n_plane_full/norm(n_plane_full);

% Defines vector xs (defined as a vector perpendicular to x_pq, within the plane, from p which intersects with the velocity vector)
x_ps = cross(n_plane,x_pq)*tan(alpha);

% Computes the velocity vector
vec_full = x_pq - x_ps;
vec = vec_full/norm(vec_full);


