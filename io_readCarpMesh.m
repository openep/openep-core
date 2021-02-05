function mesh = io_readCarpMesh(fileName,varargin)
% IO_READCARPMESH imports a Carp mesh file (.elem, .pts)
%
% Usage:
%   mesh = io_readCARPMesh(fileName)
% Where:
%   fileName - is the full path of the carp mesh
%   mesh - is a structure representing the mesh with one or more of the 
%          following fields:
%               Mesh.Pts
%                   .Tet
%                   .Hex
%                   .Tri
%                   .Quad
%                   .Edg
%                   .Pyr
%                   .Pri
%
% IO_READCARPMESH accepts the following parameter-value pairs
%   'mute'    {True} | boolean
%
% Author: Cesare Corrado (2016) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% Mesh = io_readCARPMesh(fileName);
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

file = fopen( [fileName,'.pts'], 'rt' );
if ~(mute)
    disp('Read Nodes');
end

nNod = fscanf( file, '%d', 1 );
mesh.Pts = fscanf( file, '%f %f %f', [3, nNod] )';
clear('nNod');

file = fopen( [fileName,'.elem'], 'rt' );
if ~(mute)
    disp('Read Elements');
end
Nel= fscanf( file, '%d', 1 );

%pre-allocating
%Triangles
Tri=-1.0*ones(Nel,4);
counterTri=0;

%TetraHedra
Tet=-1.0*ones(Nel,5);
counterTet=0;

%Quadrilaterals
Quad = -1.0*ones(Nel,5);
counterQuad=0;

%Hexahedra
Hex= -1.0*ones(Nel,9);
counterHex=0;

%Edges (connetions)
Edg = -1.0*ones(Nel,3);
counterEdg=0;

%pyramids
Pyr = -1.0*ones(Nel,6);
counterPyr=0;

%prisms
Pri = -1.0*ones(Nel,7);
counterPri=0;

while (~feof( file ))
    elType=fscanf(file,'%s',1);
    elType = strtrim(elType);
    switch (elType)
        case 'Tt' %Tetra
            counterTet=counterTet+1;
            data=fscanf(file,'%d %d %d %d %d',[1,5]);
            Tet(counterTet,:)=data;
            clear('data');
        case 'Tr' %Triangle
            counterTri=counterTri+1;
            data=fscanf(file,'%d %d %d %d',[1,4]);
            Tri(counterTri,:)=data;
            clear('data');            
        case 'Qd' %quadrilateral
            counterQuad=counterQuad+1;
            data=fscanf(file,'%d %d %d %d %d',[1,5]);
            Quad(counterQuad,:)=data;
            clear('data');
        case 'Hx' %Hexaedra
            counterHex=counterHex+1;
            data=fscanf(file,'%d %d %d %d %d %d %d %d %d',[1,9]);
            Hex(counterHex,:)=data;
            clear('data');
        case 'Cx' %Edges
            counterEdg=counterEdg+1;
            data=fscanf(file,'%d %d %d',[1,3]);
            Edg(counterEdg,:)=data;
            clear('data');
        case 'Py' %Pyramids
            counterPyr=counterPyr+1;
            data=fscanf(file,'%d %d %d %d %d %d',[1,6]);
            Pyr(counterPyr,:)=data;
            clear('data');
        case 'Pr' %Prisms
            counterPri=counterPri+1;
            data=fscanf(file,'%d %d %d %d %d %d %d',[1,7]);
            Pri(counterPri,:)=data;
            clear('data');
    end
end
fclose( file );

if(counterTet>0)
    Tet(:,1:4)=Tet(:,1:4)+1;
    mesh.Tet=Tet(1:counterTet,:);
end
clear('Tet','counterTet');

if(counterHex>0)
    Hex(:,1:8)=Hex(:,1:8)+1;
    mesh.Hex=Hex(1:counterHex,:);
end
clear('Hex','counterHex');

if(counterTri>0)
    Tri(:,1:3)=Tri(:,1:3)+1;
    mesh.Tri=Tri(1:counterTri,:);
end
clear('Tri','counterTri');

if(counterQuad>0)
    Quad(:,1:4)=Quad(:,1:4)+1;
    mesh.Quad=Quad(1:counterQuad,:);
end
clear('Quad','counterQuad');

if(counterEdg>0)
    Edg(:,1:2)=Edg(:,1:2)+1;
    mesh.Edg=Edg(1:counterEdg,:);
end
clear('Edg','counterEdg');

if(counterPyr>0)
    Pyr(:,1:5)=Pyr(:,1:5)+1;
    mesh.Pyr=Pyr(1:counterPyr,:);
end
clear('Pyr','counterPyr');

if(counterPri>0)
    Pri(:,1:5)=Pri(:,1:6)+1;
    mesh.Pri=Pri(1:counterPri,:);
end
clear('Pri','counterPri');


