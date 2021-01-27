function [iElectrode xyz] = read_positions_on_annotation_v2(varargin)
% READ_POSITIONS_ON_ANNOTATION_V2 loads this Carto3 position file.
% Usage:
%   [iElectrode xyz] = read_positions_on_annotation_v2(filename)
% Where:
%   iElectrode is a vector of electrode numbers
%   xyz is an array of xyz positions

% Author: Nick Linton (2013) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------


% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

filename = varargin{1};

if isempty(strfind(filename, 'OnAnnotation'))
   error('READ_POSITIONS_ON_ANNOTATION_V2: filename must be an "OnAnnotation" file.')
end

fid = fopen(filename, 'r');
if fid == (-1)
    error(['READ_POSITIONS_ON_ANNOTATION_V2: Could not read the file: "' filename '"']);
end
try
    line1 = fgetl(fid);
    line2 = fgetl(fid);
    vData = fread(fid,'*char')';
    fclose(fid);
catch err
    fclose(fid);
    rethrow(err)
end
    
%line 1
spaces = isspace(line1);
line1(spaces) = [];
line1 = lower(line1);
matchA = strfind(line1, 'eleclectrode_positions_2.0');
matchB = strfind(line1, 'sensor_positions_2.0');
if isempty(matchA) && isempty(matchB)
    warning('READ_POSITIONS_ON_ANNOTATION_V2: The version number in the txt file is unexpected.') %#ok<*WNTAG>
end

%line 2
spaces = isspace(line2);
line2(spaces) = [];
line2 = lower(line2);
matchA = strfind(line2, 'electrode#timexyz');
matchB = strfind(line2, 'sensor#timexyz');
if isempty(matchA) && isempty(matchB)
    warning('READ_POSITIONS_ON_ANNOTATION_V2: Unexpected column titles.') %#ok<*WNTAG>
end

% the rest is the data
xyz = sscanf(vData, '%f');
xyz = reshape(xyz, 5, numel(xyz)/5)';
iElectrode = xyz(:,1);
xyz = xyz(:,3:5);

