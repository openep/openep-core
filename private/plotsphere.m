function h = plotsphere(cx,cy,cz,color, r, n, varargin)
% PLOTSPHERE adds a sphere to the current axes.
% Usage:
%   h = plotsphere(cx,cy,cz,color,r,n)
%   h = plotsphere(... , 'PropertyName', 'PropertyValue' ...)
% Where:
%   h is a handle to the sphere
%   cx, cy, cz are the coordinates of the centre
%   color can be a color into the color map or true color
%   r is the radius of the sphere
%   n - see sphere(n)
% Author: Nick Linton (2006)
% Modifications - 

if not(isnumeric(n))
    error('PLOTSPHERE: n should be a positive integer')
elseif n <= 2
    n = 20;
end
[x, y, z] = sphere(n);
x = r.*x + cx;
y = r.*y + cy;
z = r.*z + cz;

color = colorspec2rgb(color);
if length(color) == 1
    c = color(ones(size(z)));
elseif length(color) == 3
    c = zeros([size(z) , 3]);
    c(:,:,1) = color(1);
    c(:,:,2) = color(2);
    c(:,:,3) = color(3);
else
    error('PLOTSPHERE: color should be a single value or an array of length 3 specifying RGB colors')
end

if nargin == 6
    h = surf(x,y,z,c);
elseif nargin > 6
    h = surf(x,y,z,c, varargin{1:end});
end
