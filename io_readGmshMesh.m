function mesh = io_readGmshMesh(fileName,varargin)
% IO_READGMSHMESH imports a GMSH mesh file (.mesh)
%
% Usage:
%   mesh = io_readGmshMesh(fileName)
% Where:
%   fileName - is the full path of the .mesh file (with the extension)
%   mesh - is a structure representing the mesh with one or more of the 
%          following fields:
%               mesh.Pts
%                   .Edg
%                   .Tri
%                   .Quad
%                   .Tet
%                   .Hex
% IO_READGMSHMESH accepts the following parameter-value pairs
%   'mute'    {True} | boolean
%
% Author: Cesare Corrado (2016) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% Mesh = io_readGmshMesh(fileName);
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 1; % UPDATE VALUE
mute=false;
if nargin > nStandardArgs
    for iArg = 1:2:nargin-nStandardArgs
        switch varargin{iArg}
            case 'mute'
                mute = varargin{iArg+1};
        end
    end
end
clear('iArg');

file = fopen( fileName, 'rt' );


while (~feof( file ))
    str = fgetl( file );
    str = strtrim(str);
    switch (str)
        case 'Vertices'
            if ~mute
                disp('Read Nodes')
            end
            nNod = fscanf( file, '%d', 1 );
            mesh.Pts = fscanf( file, '%f %f %f %f', [4, nNod] )';
        case 'Edges'
            if ~mute
                disp('Read Edges')
            end
            nEl = fscanf( file, '%d', 1 );
            mesh.Edg = fscanf( file, '%d %d %d %f', [3, nEl] )';
        case 'Triangles'
            if ~mute
                disp('Read Triangles')
            end
            nEl = fscanf( file, '%d', 1 );
            mesh.Tri = fscanf( file, '%d %d %d %f', [4, nEl] )';
        case 'Quadrilaterals'
            if ~mute
                disp('Read Quadrilaterals')
            end
            nEl = fscanf( file, '%d', 1 );
            mesh.Quad = fscanf( file, '%d %d %d %f', [5, nEl] )';
        case 'Tetrahedra'
            if ~mute
                disp('Read Tetrahedra')
            end
            nEl = fscanf( file, '%d', 1 );
            mesh.Tet = fscanf( file, '%d %d %d %d %f', [5, nEl] )';
        case 'Hexahedra'
            if ~mute
                disp('Read Hexahedra')
            end
            nEl = fscanf( file, '%d', 1 );
            mesh.Hex= fscanf( file, '%d %d %d %d %d %d %d %d %d', [9, nEl] )';
    end
end

fclose( file );
