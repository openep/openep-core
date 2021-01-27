function [l] = lineLength(h)
% LINELENGTH calculates the length of a line
%
% Usage:
%   [ l ] = lineLength( h )
% Where:
%   h  - is either
%        - a handle to a line object, or
%        - a matrix of cartesian co-ordinates representing the line, of 
%          size nx3, where n is the number of points on the line 
%
% LINELENGTH calculates the length of a line. If a handle to a line is
% passed to lineLength, then the X, Y, Z co-ordinates are taken from the
% XData, YData and ZData properties of that line object. Otherwise, the X,
% Y, Z data are received directly in a matrix of the form:
%         [ x_1  y_1  z_1 ]
%         [ x_2  y_2  z_2 ]
%         [ ...  ...  ... ]
%         [ x_n  y_n  z_n ]
%
% Author: Steven Williams (2016)
% Modifications -
%   23-09-2020: Added documentation
%   23-09-2020: Added functionality to pass in the co-ordinates directly
%
% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

%LINELENGTH Calculates the length of a line
%   h - handle to a line object
%   l - length of that line

% get the X, Y and Z values
if(ishandle(h))
    X = get(h, 'XData');
    Y = get(h, 'ZData');
    Z = get(h, 'YData');
else
    X = h(:,1);
    Y = h(:,2);
    Z = h(:,3);
end
X=X(:); Y=Y(:); Z=Z(:);

% remove any NaN values
data = [X Y Z];
data(any(isnan(data), 2), :) = [];

dX = diff(data(:,1));
dY = diff(data(:,2));
dZ = diff(data(:,3));

l = 0;
for i = 1:length(dX)
    l = l + sqrt(dX(i)^2 + dY(i)^2 + dZ(i)^2);
end

end