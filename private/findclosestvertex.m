function [vertices, distances] = findclosestvertex(t, loc, varargin)
% FINDCLOSESTVERTEX searches to find the closest vertex to a number of locations.
% Usage:
%   [vertices, distances] = findclosestvertex(t, loc, limitToTriangulation)
%   [vertices, distances] = findclosestvertex(x, loc)
% Where:
%   t is a triangulation or TriRep object
%   x is a n*3 matrix where each row contains the coordinates of a vertex.
%   loc is [x, y, z] - the coordinates of the points we're interested in
%   limitToTriangulation - if there are points in t.Points that are not in the
%       triangulation, then set limitToTriangulation = true. This will
%       create a temporary copy of the triangulation and check it (making
%       it quite slow).
% Author: Nick Linton (2009)
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------
                        % if true
                        %     nPoints = 20000;
                        %     p = rand(nPoints,3)*2*pi;
                        %     epicentre = round(rand(1)*(nPoints-1) + 1);
                        %     for i = 1:nPoints
                        %         p(i,:) = 100* [1 1 1] * rotationmatrix(p(i,1), p(i,2), p(i,3));
                        %     end
                        %     tri = convhulln(p);
                        %     t = TriRep(tri,p);
                        %     testpoints = rand(5,3);
                        % end
                        % profile on
                        % profile clear
                        % for i = 1:5
                        %     vertex = findclosestvertex(t, testpoints(i,:))
                        % end
                        % profile viewer

% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

    % handle the limitToTriangulation
    if nargin == 3 && (isa(t,'TriRep') || isa(t,'triangulation'))
        if varargin{1}==true
            warning('FINDCLOSESTVERTEX:hasToTrim' , 'FINDCLOSESTVERTEX: runs slowly if it has to trim a triangulation') 
            [tNew, isVertUsed] = repack(t);
            [verticesNew, distances] = findclosestvertex(tNew, loc);
            iNewVertex = cumsum(isVertUsed);
            vertices = zeros(size(verticesNew));
            for i = 1:numel(verticesNew)
                vertices(i) = find(iNewVertex==verticesNew(i),1,'first');
            end
            return
        end
    end
    
    if isa(t, 'double')
        % convert it to a structure, similar to a triangulation to fool the code below.
        temp.Points = t;
        t = temp;
        if nargin~=2
            error('FINDCLOSESTVERTEX: wrong number of arguments for vertices input (rather than TriRep)')
        end
    elseif isa(t,'TriRep')
        temp.Points = t.X;
        t = temp;
    elseif isa(t,'triangulation')
        %then nothing to do
    else
        error('FINDCLOSESTVERTEX: wrong class of input for first argument.')
    end

    vertices = zeros(size(loc,1),1);
    distances = zeros(size(vertices));
    dist = 0;
    for i = 1:size(loc,1)
        loc2 = loc(i*ones(size(t.Points, 1),1), :); % similar to repmat but quicker as less checking
        test = (t.Points - loc2);
        test = sum(test.*test, 2);
        if nargout == 1
            [~,vertex] = min(test); %this is the fastest version
        elseif nargout == 2
            [dist,vertex] = min(test);
            dist = sqrt(dist);
        else
            error('FINDCLOSESTVERTEX: Wrong number of output arguments')
        end
        vertices(i) = vertex;
        distances(i) = dist;
    end
    
end