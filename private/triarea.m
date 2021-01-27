function area = triarea(varargin)
% TRIAREA calculates the areas in a triangulation.
% Usage:
%   area = triarea(trirep)
%   area = triarea(simplices, points)
%   area = triarea(points)
% Where:
%   area - area of each triangle in trirep, triangulation, or a single value
%   trirep - a TriRep object
%   triangulation - a triangulation object
%   simplices and points - the parameters as per format of TriRep
%   points - a 3*3 or 3*2 matrix of points comprising a single triangle
% TRIAREA computes the area of the triangles using Heron's formula.

% Author: Nick Linton (2010)
% Modifications - 2013 capability for triangulation object added

switch nargin
    case 1
        if isa(varargin{1},'TriRep')
            simplices = varargin{1}.Triangulation;
            points = varargin{1}.X;
        elseif isa(varargin{1},'triangulation')
            simplices = varargin{1}.ConnectivityList;
            points = varargin{1}.Points;
        elseif isa(varargin{1},'double')
            simplices = [1, 2, 3];
            points = varargin{1};
        else
            error('Unexpected format of input variables')
        end
    case 2
        simplices = varargin{1};
        points = varargin{2};
end

if not(size(points,2)==3 || size(points,2)==2)
    error('"points" should be a n*2 or n*3 matrix')
end
if size(simplices,2)~=3
    error('"triangulation" should be a n*3 matrix')
end

% The version of Heron's formula used only uses one sqrt for speed.
a = points(simplices(:,1),:) - points(simplices(:,2),:);
b = points(simplices(:,1),:) - points(simplices(:,3),:);
c = points(simplices(:,3),:) - points(simplices(:,2),:);

a2 = sum(a.*a,2); % (length of a)^2
b2 = sum(b.*b,2);
c2 = sum(c.*c,2);

s2 = a2+b2+c2;                  % s2 = sum of squared lenths
s4 = a2.*a2 + b2.*b2 + c2.*c2;  % s4 = sum of fourth powered lengths

area = 0.25 * sqrt(  s2.*s2 - 2*s4  );
