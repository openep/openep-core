function D = distBetweenPoints(A, B, varargin)
%     **********************************************************************
%     * The contents of this package are copyright (c) King's College London
%     *
%     * They may not be reproduced, distributed, modified or sold for any 
%     * purpose
%     *
%     * Author: 	Dr S. E. Williams
%     * Address: 	Division of Imaging Sciences and Biomedical Engineering,
%     *          	King's College London
%     *			St Thomas' Hospital
%     *
%     * March 2014
%     **********************************************************************
%
% DISTBETWEENPOINTS Returns the distance from A to B.
% Usage:
%   D = DISTBETWEENPOINTS(A, B)
% Where:
%   A - is the first point(s)
%   B - is the second point(s)
%
% DISTBETWEENPOINTS returns the distance from A to B. A and B are specified
% as row vectors [x, y, z] or matrices, with rows representing different
% points. If npoints in A and B are different A must specify one and only 
% one point.
%
% Author: Steven Williams (2013)
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

if nargin==3
    warn = varargin{1};
else
    warn = false;
end
if size(A,2) ~= size(B,2);
    error('DISTBETWEENPOINTS: dimensions of points in A and B must be equal');
end
if size(A,1) ~= size(B,1);
    if size(A,1) ~= 1
        error('DISTBETWEENPOINTS: if npoints in A and B are different A must specify one and only one point');
    else
        A = repmat(A, size(B,1), 1);
        if warn
            warning('DISTBETWEENPOINTS: A has been replicated size(B,1) times');
        end
    end
end

diffsq = (A - B).^2;

D = sqrt(sum(diffsq, 2));

end
