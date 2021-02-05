function writeOpenEP2CarpMesh(userdata, rootMeshFile, varargin)
% WRITEOPENEP2CARPMESH exports a mesh structure to a Carp format (.elem, .pts)
%
% Usage:
%   writeOpenEP2CarpMesh(userdata, rootMeshFile)
% Where:
%   userdata - see importcarto_mem
%   rootMeshFile - is the name (suffix) of the output carp mesh
%
% WRITEOPENEP2CARPMESH accepts the following parameter-value pairs
%   'mute'    {True} | boolean
%
% Author: Steven Williams (2021) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% writeOpenEP2CarpMesh(userdata, rootMeshFile);
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 2; % UPDATE VALUE
mute=false;
if nargin > nStandardArgs
    for iArg = 1:2:nargin-nStandardArgs
        switch varargin{iArg}
            case 'mute'
                mute = varargin{iArg+1};
        end
    end
end

% create a mesh structure from openep data
mesh.Pts = userdata.surface.triRep.X;
mesh.Elem = userdadta.surface.triRep.Triangulation;

% call io_writeCarpMesh
io_writeCarpMesh(rootMeshFile, mesh, varargin);