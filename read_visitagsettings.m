function [visitagSettings] = read_visitagsettings(varargin)
% READ_VISITAGSETTINGS loads Carto3 visitag settings file.
% Usage:
%   [visitagSettings] = read_visitagsettings(filename)
% Where:
%   visitagSettings is a structure of the visitag settings
%   filename is the full file path
%
% READ_VISITAGSETTINGS creates a structure with the parameters in
% VisiTagSettings.txt. Data is converted to double if it is numeric or
% remains as a string if it is a string. Beware of equals sign at the end
% of the parameter names - these are currently removed by the code but
% would be an obvious source of future errors.
%
% Author: Steven Williams (2015) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------
% visitagSettings = read_visitagsettings(filename);
% ---------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

filename = varargin{1};

if isempty(strfind(filename, 'VisiTagSettings'))
    error('READ_VISITAGSETTINGS: filename must be a "VisiTagSettings" file.')
end

fid = fopen(filename, 'r');
if fid == (-1)
    error(['READ_VISITAGSETTINGS: Could not read the file: "' filename '"']);
end
try
    tLine = fgetl(fid);
    while ischar(tLine)
        C = strsplit(tLine);
        C(strcmpi('',C)) = [];
        param = str2double(C{2});
        if isnan(param)
            param = C{2};
        end
        visitagSettings.(C{1}(1:end-1)) = param;
        tLine = fgetl(fid);
    end
    fclose(fid);
catch
    fclose(fid);
    rethrow(err);
end

end