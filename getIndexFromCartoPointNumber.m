function [index] = getIndexFromCartoPointNumber(userdata, pointNumber)
% GETINDEXFROMCARTOPOINTNUMBER Finds the index of the mapping point at the 
% point number displayed on the Carto mapping system.
%
% Usage:
%   [index] = getIndexFromCartoPointNumber(pointNumber)
% Where:
%   userdata   - a userdata structure
%   pointNumber - a point number (or array of point numbers) as displayed
%                 on the Carto mapping system
%   index       - an index (or array of indices) for referencing into the
%                 data fields within userdata.electric
%
% GETINDEXFROMCARTOPOINTNUMBER Finds the index of the mapping point at the 
% point number displayed on the Carto mapping system.
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications - 
%
% Info on Code Testing:
% ---------------------------------------------------------------
% index = getIndexFromCartoPointNumber(userdata, 1)
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

index = [];
for i = 1:numel(pointNumber)
    point = ['P' num2str(pointNumber(i))];
    T = find(strcmpi(userdata.electric.names, point));
    if isempty(T)
        index(i) = NaN;
        warning(['GETINDEXFROMCARTOPOINTNUMBER: Point ' point ' does not exist'])
    else
        index(i) = T;
    end
end

end



          

    
    