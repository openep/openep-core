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
            for i = 1:nVertices
                tLine = fgetl(fid);
                endText = find(tLine=='=',1);
                temp = str2num(tLine((endText+1):end)); %#ok<ST2NM>
                xyz(i,:) = temp(1:3);
                normals(i,:) = temp(4:6);
            end
            %check that the next line is blank
            tLine = fgetl(fid);
            if ~isempty(tLine); error('READ_MESHFILE: unexpected format.'); end
        elseif strstartcmpi('[TrianglesSection]', tLine)
            % there is one line of headers then a blank line
            tLine = fgetl(fid); %#ok<NASGU>
            tLine = fgetl(fid);
            if ~isempty(tLine); error('READ_MESHFILE: unexpected format.'); end
            
            tri = zeros(nTriangles, 3);
            groupId = zeros(nTriangles,1);
            for i = 1:nTriangles
                tLine = fgetl(fid);
                endText = find(tLine=='=',1);
                temp = str2num(tLine((endText+1):end)); %#ok<ST2NM>
                tri(i,:) = temp(1:3);
                groupId(i) = temp(7);
            end
            tri = 1+tri; %convert to 1-based indexing
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
            
            act_bip = zeros(nVertices, 2);
            uni_imp_frc = zeros(nVertices, 3);
            for i = 1:nVertices
                tLine = fgetl(fid);
                endText = find(tLine=='=',1);
                temp = str2num(tLine((endText+1):end)); %#ok<ST2NM>
                temp(temp == -10000) = NaN;
                act_bip(i,:) = temp([3,2]);
                if nargout>=5
                    uni_imp_frc(i,:) = temp([1,4,11]);
                end
            end

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