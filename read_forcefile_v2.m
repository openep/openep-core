function varargout = read_forcefile_v2(varargin)
% READ_FORCEFILE loads this Carto3 force file.
% Usage:
%   force = read_forcefile(filename);
%   [force axialAngle lateralAngle] = read_forcefile(filename)
%   [force axialAngle lateralAngle t_time t_force t_axialAngle t_lateralAngle] = read_forcefile(filename)
% Where:
%   force - is the fixed time point force for this point
%   axialAngle - is the fixed time point axial angle for this point
%   lateralAngle - is the fixed time point lateral angle for this point
%   t_time - is the time array (e.g. -7000ms->5000ms) for the time data
%   t_force - is the time course of force
%   t_axialAngle - is the time course of the axial angle
%   t_lateralAngle - is the time course of the lateral angle
%
%   filename is the Carto3.txt file
%
% Author: Steven Williams (2013) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications - 
%
% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------

% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% load data from file
filename = varargin{1};
fid = fopen(filename, 'r');
if fid == (-1)
    error(['READ_FORCEFILE: Could not read the file: "' filename '"']);
end
try
    line1 = fgetl(fid);
    line2 = fgetl(fid);
    line3 = fgetl(fid);
    line4 = fgetl(fid); % Tyically IntervalGraph=100
    line5 = fgetl(fid); % Contains: Time, Force, AxialAngle, LateralAngle, MetalSeverity 
    line6 = fgetl(fid); % Typically IntervalNonGraph=1000
    line7 = fgetl(fid); % Contains: Time, Force, AxialAngle, LateralAngle, MetalSeverity
    line8 = fgetl(fid);
    if ~strcmpi(line4(1:13),'IntervalGraph')
            error(['READ_FORCEFILE: Could not read the file: "' filename '"'])
    end
    if ~strcmpi(line8(1:5),'index')
        %then we may have another version of the ecgfile
        %try the next line
        line8 = fgetl(fid);
        if ~strcmpi(line8(1:5),'index')
            error(['READ_FORCEFILE: Could not read the file: "' filename '"'])
        end
    end
    fData = fread(fid,'*char')';
    fclose(fid);
catch err
    fclose(fid);
    rethrow(err)
end
    
%line 1
spaces = isspace(line1);
line1(spaces) = [];
line1 = lower(line1);
match = strfind(line1, 'contactforce.txt_2.0');
if isempty(match)
    warning('READ_FORCEFILE: The version number in the txt file is unexpected.') %#ok<*WNTAG>
end

%line 2
spaces = isspace(line2);
line2(spaces) = [];
line2 = lower(line2);
match = regexp(line2, '\d+', 'match');
numTimePoints = str2double(match{2});
if ~strcmpi(match{1},'50') || ~strcmpi(match{2},'200')
    warning(['READ_FORCEFILE: Unexpected rate (found: ' match{1} ' expected: 50) or number (found: ' match{2} ' expected: 200)']) %#ok<*WNTAG>
end

%line 3
spaces = isspace(line3);
line3(spaces) = [];
line3 = lower(line3);
match = strfind(line3, 'mode=0');
if isempty(match)
    warning(['READ_FORCEFILE: Unexpected mode found']);
end

%line 4
% ignore line 4

%line 5
match = regexp(line7, '\d+\.?\d*', 'match');
force = match{2};
axialAngle = match{3};
lateralAngle = match{4};

if nargout == 1
    varargout{1} = force;
elseif nargout == 3
    varargout{1} = force;
    varargout{2} = axialAngle;
    varargout{3} = lateralAngle;
elseif nargout == 8 % the user has asked for timecourse data
    varargout{1} = force;
    varargout{2} = axialAngle;
    varargout{3} = lateralAngle;
    % deal with the timecourse data
    fD = reshape(sscanf(fData, '%f'), [9 numTimePoints])';
    % the 9 columns of fD are in column order:
    %   INDEX TIME TIMESTAMP FORCE AXIALANGLE LATERALANGLE METALSEVERITY
    %   INACCURATEMETALSEVERITY NEEDZEROING
    varargout{4} = fD(:,2); % TIME
    varargout{5} = fD(:,4); % FORCE
    varargout{6} = fD(:,5); % AXIALANGLE
    varargout{7} = fD(:,6); % LATERALANGLE
    varargout{8} = fD(:,3); % systemtime
else
    error('READ_FORCEFILE: wrong number of output arguments.')
end

end

