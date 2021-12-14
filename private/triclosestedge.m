function [loc, dist, edgeIndex, frac] = triclosestedge(point, tri)
% TRICLOSESTEDGE finds the closest edge on a triangulation.
% Usage:
%   [loc, dist, edgeIndex, frac] = triclosestedge(point, tri)
% Where:
%   loc is the closest location on the triangulation edge
%   dist is distance from point to triangulation edge
%   edgeIndex gives the edge row in edges(tri)
%   frac gives the distance along edge(1) edge(2) as a proportion (0 to 1)
%   point is a point [x,y,z], or can be an array of points.
%   tri is a triangulation
%
% TRICLOSESTEDGE finds the closest location on the triangulation edges.
%
% Author: Nick Linton (2017)

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------

                        

% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------
   
    % get the triangulation info
    if isa(tri, 'triangulation')
        x = tri.Points;
        cl = tri.ConnectivityList;
        ed = edges(tri);
    elseif isa(tri, 'TriRep')
        x = tri.X;
        cl = tri.Triangulation;
        ed = edges(tri);
    else
        error('Wrong input type')
    end
    
    
    
    % sort out the data for point
    validateattributes(point, {'numeric'}, {'ndims',2});
    if ~isvector(point)
        if size(point,2)~=size(x,2)
            error('TRICLOSESTEDGE: point has different dimensions to the triangulation.');
        end
        n = size(point,1);
        loc = zeros(size(point));
        dist = zeros(n,1);
        frac = zeros(n,1);
        edgeIndex  = zeros(n,1);
        for i = 1:n
            [loc(i,:), dist(i), edgeIndex(i), frac(i)] = local_triclosestedge(point(i,:), x, ed);
        end
    else
        [loc, dist, edgeIndex, frac] = local_triclosestedge(point(:)', x, ed);
    end
    
end

function [loc, dist, iEdge, frac] = local_triclosestedge(point, x, ed)
    
    % To avoid multiple subtractions, translate everything so that point is
    % at the origin.
    x = bsxfun(@minus, x, point);
    
    % start with a first guess
    distSq = sum(x.*x,2);
    [bestDistSq, v] = min(distSq);  % to avoid too many sqrt operations;

    if bestDistSq == 0
        loc = x(v,:) + point;
        % take first edge that has v as one of the vertices
        for iEd = 1:size(ed,1)
            if ed(iEd,1) == v
                iEdge = iEd; dist = 0; frac = 0; return
            elseif ed(iEd,2) == v
                iEdge = iEd; dist = 0; frac = 1; return
            end
        end
    else
        % We can immediately reject all edges for which both points are more
        % than dist along a particular direction - for example the direction
        % from point to the middle of the edge in question.
        
        % To avoid a division calculation (for speed), let p1 and p2 be
        % points on an edge. pM is the vector to their midpoint. nM is the
        % normalised vector to their midpoint.
        % if p1.nM > bestDist  and p2.nM > bestDist then exclude this edge
        % or
        % if (p1.pM)^2 > bestDistSq * |pM|^2   then exclude this edge
        
        isPossibleEdge = true(size(ed,1),1);
        edIndices = 1:size(ed,1);
        pM = 0.5 * (x(ed(:,1),:) + x(ed(:,2),:));
        pM_MagSq = sum(pM.*pM,2);
        
        % taking the first vertex of each simplex
        p1_dot_pM = sum( x(ed(:,1),:) .* pM  ,  2);
        % taking the second vertex of each simplex
        p2_dot_pM = sum( x(ed(:,2),:) .* pM  ,  2);
        
        temp = bestDistSq * pM_MagSq;
        
        result1 =  (p1_dot_pM.*p1_dot_pM > temp);
        result2 =  (p2_dot_pM.*p2_dot_pM > temp);
       % isPossibleEdge(result1 & result2) = false;

        ed = ed(isPossibleEdge,:);    edIndices = edIndices(isPossibleEdge(:));
        
        % Test to get all the closest points to the faces:
        [newDistSq, closestPoint, frac] = local_origin2line( x(ed(:,1),:), x(ed(:,2),:) );
        [newDistSq, iMin] = min(newDistSq);

        bestDistSq = newDistSq;
        loc = closestPoint(iMin,:) + point;
        iEdge = edIndices(iMin);
        frac = frac(iMin);

        
    end
    
    dist = sqrt(bestDistSq);

end
            
function [distSq, closestPoint, x] = local_origin2line( lineA, lineB )
%                    origin
%                     /|
%                  b / |
%                   /  |d
%               lineA------------------lineB
%                         a
% d = -b + x*a%x unknown, 0<x<1 if d ends up between lin1 and line2
% dot(a,d) = 0
% dot(a,(-b+x*a)) = 0  -> x = dot(b,a)/dot(a,a)

    a = lineB - lineA;
    b = - lineA; %origin - lineA
    x = sum(b.*a,2) ./ sum(a.*a,2);
    
    %restrict to between points:
    x(x<0) = 0;
    x(x>1) = 1;

    closestPoint = lineA + x(:,ones(1,size(lineA,2))) .* a;
    distSq = sum(closestPoint.*closestPoint,2);
end
            