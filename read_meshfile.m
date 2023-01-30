function varargout = read_meshfile(filename)
% READ_MESHFILE loads this Carto3 mesh file.
% Usage:
%   t = read_meshfile(filename)
%   [t isVertexAtEdge] = read_meshfile(filename)
%   [t isVertexAtEdge act_bip normals] = read_meshfile(filename)
% Where:
%   filename is the filename
%   t is a TriRep object
%   isVertexAtEdge is a logical array indicating vertices at an edge
%   act_bip is a matrix of activation times and bipolar voltage by vertex
%
% READ_MESHFILE reads the triangulation from a TriangulatedMeshVersion 2.0
% file.

% Author: Nick Linton (2011) (Copyright)
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

t = [];
isVertexAtEdge = [];
act_bip = [];
normals = [];
uni_imp_frc = [];

fid = fopen(filename, 'r');
if fid == (-1)
    warning(['READ_MESHFILE: Could not open file: ' filename]);
    beep()
    varargout = cell(1,nargout);
    return
end
cleanupObj = onCleanup(@()fclose(fid));

    %line 1
    tLine = fgetl(fid);
    spaces = isspace(tLine);
    tLine(spaces) = [];
    tLine = lower(tLine);
    match = strfind(tLine, 'triangulatedmeshversion2.0');
    if isempty(match)
        warning('READ_MESHFILE: The version number in the txt file is unexpected.') %#ok<WNTAG>
    end
    eof = false;

    while (~eof)
        tLine = fgetl(fid);
        if tLine == (-1)
            eof = true;
        elseif strstartcmpi('numvertex ', tLine)
            endText = find(tLine=='=',1);
            nVertices = str2double(tLine((endText+1):end));
        elseif strstartcmpi('numtriangle ', tLine)
            endText = find(tLine=='=',1);
            nTriangles = str2double(tLine((endText+1):end));
        elseif strstartcmpi('[VerticesSection]', tLine)
            % there is one line of headers then a blank line
            tLine = fgetl(fid); %#ok<NASGU>
            tLine = fgetl(fid);
            if ~isempty(tLine); error('READ_MESHFILE: unexpected format.'); end
            xyz = zeros(nVertices, 3);
            normals = zeros(nVertices, 3);
            
            formatSpec = '%d = %f %f %f %f %f %f %d';
            data = fscanf(fid,formatSpec,[nVertices*8 , 1]);
            data = reshape(data, 8,nVertices)';
            
            xyz = data(:,2:4);
            normals = data(:,5:7);
            groupId = data(8);
            
            %check that the next line is blank
            tLine = fgetl(fid);
            if ~isempty(tLine); error('READ_MESHFILE: unexpected format.'); end
        elseif strstartcmpi('[TrianglesSection]', tLine)
            % there is one line of headers then a blank line
            tLine = fgetl(fid); %#ok<NASGU>
            tLine = fgetl(fid);
            if ~isempty(tLine); error('READ_MESHFILE: unexpected format.'); end
            
            formatSpec = '%d = %d %d %d %f %f %f %d';
            data = fscanf(fid,formatSpec,[nTriangles*8 , 1]);
            data = reshape(data, 8,nTriangles)';
            
            tri = data(:,2:4)+1; %1-based indexing
            groupId = data(8);
            %normalVector = data(:,5:8);
            tri(groupId<0,:) = [];
            
            warning('off', 'MATLAB:TriRep:PtsNotInTriWarnId')
                t = TriRep(tri, xyz);
            warning('on', 'MATLAB:TriRep:PtsNotInTriWarnId')
            
            % Now find the vertices which make up the edges. We find the
            % faces which are at the freeBoundary. Taking all of the
            % vertices, a vertex which is at the edge must be listed twice
            % (the way that the Carto anatomy is cut out).
            isVertexAtEdge = false(size(xyz,1),1);
            vFree = freeBoundary(t);
            vFree = vFree(:);
            [~, iFB, ~] = unique(vFree);
            vFree(iFB) = [];    % we have wiped out all vertices which are only listed once
            isVertexAtEdge(vFree) = true;

            %check that the next line is blank
            tLine = fgetl(fid);
            if ~isempty(tLine); error('READ_MESHFILE: unexpected format.'); end
            
        elseif strstartcmpi('[VerticesColorsSection]', tLine)
            % there are 2 lines of headers then a blank line
            tLine = fgetl(fid); %#ok<NASGU>
            tLine = fgetl(fid); %#ok<NASGU>
            tLine = fgetl(fid);
            if ~isempty(tLine); error('READ_MESHFILE: unexpected format.'); end
            
            formatSpec = '%d = %f %f %f %f %f %f %f %f %f %f %f %f %f';
            data = fscanf(fid,formatSpec,[nVertices*14 , 1]);
            data = reshape(data, 14,nVertices)';
            
            data(data==-10000) = NaN;
            act_bip = data(:,[4,3]);
            uni_imp_frc = data(:,[2,5,12]);
            
            %check that the next line is blank
            tLine = fgetl(fid);
            if ~isempty(tLine); error('READ_MESHFILE: unexpected format.'); end
            break
        end
    end
    
varargout{1} = t;
if nargout>=2; varargout{2} = isVertexAtEdge;   end
if nargout>=3; varargout{3} = act_bip;          end
if nargout>=4; varargout{4} = normals;          end
if nargout>=5; varargout{5} = uni_imp_frc;      end