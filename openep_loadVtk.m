function userdata = openep_loadVtk( varargin )
% OPENEP_LOADVTK Loads a VTK file and creates an OpenEP dataset
%
% Usage:
%   userdata = openep_loadVtk( varargin )
% Where:
%   varargin - a path to a VTK file; or 'openfile')
%   userdata  - an OpenEP dataset
%
% OPENEP_LOADVTK accepts the following parameter-value pairs
%   'param1'     {value1}|vallue2
%
% OPENEP_LOADVTK Detailed description goes here
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

in1 = varargin{1};
if strcmpi(in1, 'openfile')
    formatspec = {'*.vtk'};
    instructionstring = 'Select a file';
    [filename, pathname] = uigetfile(formatspec, instructionstring);
    filePath = [pathname filename];
else
    filePath = in1;
end

hVtk = VTKReader(filePath);
hVtk.readAllData();

userdata = openep_createuserdata();
tr = TriRep(hVtk.Cells+1, hVtk.Points); % +1 for indexing

userdata.surface.triRep = tr;

end