function out = convertDataToMap(userdata, X, data)
% CONVERTDATATOMAP maps data at known co-ordinates to the surface map
%
% Usage:
%   out = convertDataToMap(userdata, X, data)
% Where:
%   userdata - an OpenEP data structure
%   X        - locations of data, size n x 3
%   data     - scalar values, size n x 1
%
% CONVERTDATATOMAP accepts the following parameter-value pairs
%   'param1'     {value1}|vallue2
%
% CONVERTDATATOMAP Detailed description goes here
%
% Author: Steven Williams (2022)
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% identify vertices, quantify distances
[iV, distances] = findclosestvertex(getMesh(userdata), X, true);
meanDistance = mean(distances);
stdDistances = std(distances);

% output to command window
disp(['Found locations for ' ...
    num2str(numel(iV)) ...
    ' vertices. Averate distance ' ...
    num2str(meanDistance) ...
    '±' ...
    num2str(stdDistances) ...
    '.' ...
    ])

% assign output
out = NaN(length(getVertices(userdata)),1);
out(iV) = data;

end