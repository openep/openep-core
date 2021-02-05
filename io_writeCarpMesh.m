function io_writeCarpMesh(rootMeshFile,mesh,varargin)
% IO_WRITECARPMESH exports a mesh structure to a Carp format (.elem, .pts)
%
% Usage:
%   io_writeCarpMesh(rootMeshFile,mesh)
% Where:
%   rootMeshFile - is the name (suffix) of the output carp mesh
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
% IO_WRITECARPMESH accepts the following parameter-value pairs
%   'mute'    {True} | boolean
%
% Author: Cesare Corrado (2016) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% io_writeCARPMesh(rootMeshFile, mesh);
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
clear('iArg');

meshFields{1}='Tet';
meshFields{2}='Hex';
meshFields{3}='Tri';
meshFields{4}='Quad';
meshFields{5}='Pyr';
meshFields{6}='Pri';
meshFields{7}='Edg';

ElemName{1}='Tt';
ElemName{2}='Hx';
ElemName{3}='Tr';
ElemName{4}='Qd';
ElemName{5}='Py';
ElemName{6}='Pr';
ElemName{7}='Cx';

ElemNbNodes(1)=4;
ElemNbNodes(2)=8;
ElemNbNodes(3)=3;
ElemNbNodes(4)=4;
ElemNbNodes(5)=5;
ElemNbNodes(6)=6;
ElemNbNodes(7)=2;

%write nodes
fidend = fopen([rootMeshFile,'.pts'],'w');
nPt=size(mesh.Pts,1);
if ~mute
    disp('Write Points')
end
fprintf(fidend,'%d',nPt);
for iPt=1:nPt
    fprintf(fidend,'\n%-10.6f %-10.6f %-10.6f',mesh.Pts(iPt,1:3));
end
fclose(fidend);
clear('fidend','iPt','nPt');

nElem=0;
shift=true;
for iField=1:numel(meshFields)
    if(isfield(mesh,meshFields{iField}))
        nElem=nElem+size(mesh.(meshFields{iField})  ,1);
        Elements=mesh.(meshFields{iField});
        Elements=Elements(:,1:ElemNbNodes(iField)  );
        if(int32(min( Elements(:) )  )==0)
            shift=false;
        end
    end
end

%write Elements
fidend = fopen([rootMeshFile,'.elem'],'w');
if ~mute
    disp('Write Elements')
end
fprintf(fidend,'%d',nElem);
clear('nElem');
for iField=1:numel(meshFields)
    if(isfield(mesh,meshFields{iField}))
        Elements=mesh.(meshFields{iField});
        [nElem,nEntry]=size(Elements);
        if(shift)
            Elements(:,1:ElemNbNodes(iField))=Elements(:,1:ElemNbNodes(iField))-1;
        end
        for iElem=1:nElem
            fprintf(fidend,'\n%s',ElemName{iField});
            for iEntry=1:nEntry
                fprintf(fidend,' %d',Elements(iElem,iEntry));
            end
        end
        clear('Elements','nElem','nEntry','iElem','iEntry');
    end %end  on if isfield
end  %end for on fields
fclose(fidend);
clear('fidend','jField');

