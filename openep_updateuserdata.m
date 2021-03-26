function path = openep_updateuserdata(varargin)
% OPENEP_UPDATEUSERDATA is used to make non-data related changes to
% userdata
%
% Usage:
%   path = openep_updateuserdata(userdata)
%   path = openep_updateuserdata(dir)
%   path = openep_updateuserdata(f)
%
% Where:
%   userdata - openEP userdata structure
%   dir      - path to a userdata .mat file or a directory containing multiple
%              userdata .mat files
%   path     - path to the new .mat file
%
% OPENEP_UPDATEUSERDATA accepts the following parameter-value pairs
%   'version'     {'-v7.3'} | '-v7' | '-v6'
%
% OPENEP_UPDATEUSERDATA is used to make non-data related changes to
% userdata. The following interventions are possible
% 1) By specifying 'version', the save version can be set (see help save)
%
% Author: Steven Williams (2021)
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

nStandardArgs = 2; % UPDATE VALUE
value1 = 'NaN';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'param1'
                value1 = varargin{i+1};
        end
    end
end

if strcmpi(in1, 'openfile')
    formatspec = {'*.vtk'};
    instructionstring = 'Select a file';
    [filename, pathname] = uiputfile(formatspec, instructionstring);
    filePath = [pathname filename];
else
    filePath = in1;
end



end