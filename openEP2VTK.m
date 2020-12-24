function path2VTKfile = openEP2VTK(userdata, varargin)
% OPENEP2VTK Converts an OpenEP data structure to a VTK file
%
% Usage:
%   tr = openEP2VTK('openfile')
% Where:
%   tr  - a TriRep object
%   path2VTKfile - the path to the file that was written
%
% OPENEP2VTK accepts the following parameter-value pairs
%   'datatype'     {'bip'} | 'uni' | 'lat'
%       - the required data, bipolar voltage, unipolar voltage or local
%       activation time
%   'method'       {'map'} | 'egm'
%       - the method of accessing the data; clinical-system map based
%       ('map') or re-inteprolated by OpenEP from the raw egms ('egm');
%   'outputfile'
%       - absolute path to the output file. 
%       If empty then openEP2VTK outputs the VTK file to the current 
%        directory using the current date/time as the filename
%       If strcmpi('outputfile','openfile') a GUI is used to place the
%        file
%
% OPENEP2VTK Converts between OpenEP format and VTK format. This
% function takes map data and writes it to the VTK file, or if 'method' is
% set to 'egm' it first uses generateInterpData.m to create interpolated
% data.
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% path2VTKfile = openEP2VTK(userdata, 'datatype', 'lat', 'outputfile', 'openfile')
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% parse input arguments
nStandardArgs = 1;
datatype = 'lat';
method = 'map';
outputfile = '';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'datatype'
                datatype = varargin{i+1};
            case 'method'
                method = varargin{i+1};
            case 'outputfile'
                method = varargin{i+1};
        end
    end
end

% get the file location
if isempty(outputfile)
    path2VTKfile = datestr(datetime());
    path2VTKfile(isspace(path2VTKfile)) = '_';
    path2VTKfile(path2VTKfile==':') = '_';
    path2VTKfile = [path2VTKfile '.vtk'];
elseif strcmpi(outputfile, 'openfile')
    formatspec = {'*.vtk'};
    instructionstring = 'Select a file';
    [path2VTKfile, pathname] = uiputfile(formatspec, instructionstring);
    path2VTKfile = [pathname path2VTKfile];
else
    path2VTKfile = outputfile;
end

switch lower(datatype)
    case 'bip'
        switch lower(method)
            case 'map'
                pointdata = userdata.surface.act_bip(:,2);
            case 'egm'
                disp('Interpolating data ...');
                pointdata = generateInterpData(userdata, 'bip-map');
        end
        
    case 'uni'
        switch lower(method)
            case 'map'
                pointdata = userdata.surface.uni_imp_frc(:,1);
            case 'egm'
                disp('Interpolating data ...');
                pointdata = generateInterpData(userdata, 'uni-map');
        end
    case 'lat'
        switch lower(method)
            case 'map'
                pointdata = userdata.surface.act_bip(:,1);
            case 'egm'
                disp('Interpolating data ...');
                pointdata = generateInterpData(userdata, 'uni-map');
        end
end

writeTriRep2VTK(getMesh(userdata), pointdata, 'outputfile', path2VTKfile);

end
