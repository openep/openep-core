function [tNew, isVertUsed] = repack(t)
%REPACK deletes vertices that aren't referenced by a triangulation/TriRep.
% Usage:
%   tNew = repack(t)
%   [tNew, isVertUsed] = = repack(t)
% Where:
%   t is a TriRep
%   iVert is logical and t.Points(isVertUsed)=tNew.Points
%
% REPACK discards unused vertices. Note that the vertices that are kept
% will change position in the X matrix.
%
% Author: Nick Linton (2011)
% Modifications - 
%   NL 2013: altered to cope with Triangulation and to return isVertUsed.

    if isa(t,'TriRep')
        tri = t.Triangulation;
        x = t.X;
    elseif isa(t,'triangulation')
        tri = t.ConnectivityList;
        x = t.Points;
    elseif isa(t, 'struct')
        tri = t.Triangulation;
        x = t.X;
    else
        error('REPACK: wrong type of input.')
    end

    nV = size(x,1);

    isVertUsed = false(nV,1);
    isVertUsed(tri(:)) = true;

    newX = x(isVertUsed,:);
    newIndex = cumsum(double(isVertUsed));

    newTri = newIndex(tri);
    
    
    if isa(t,'TriRep')
        tNew = TriRep(newTri, newX); %#ok<DTRIREP>
    elseif isa(t,'triangulation')
        tNew = triangulation(newTri, newX);
    elseif isa(t, 'struct')
        tNew.X = newX;
        tNew.Triangulation = newTri;
    end


end



