function tr = getOpenEPSubsetFromImgMap( userdata, mOpenEP3D, l, varargin)
% GETOPENEPSUBSETFROMIMGMAP Splits the triangulation in userdata based on
% the specified label, l, referencing into mOpenEP3D
%
% Usage:
%   tr = getOpenEPSubsetFromImgMap( userdata, mOpenEP3D, l)
% Where:
%   userdata  - see importcarto_mem
%   mOpenEP3D - a mesh structure with .Pts and .Tri fields, where the
%               fourth column of .Tri is the labels
%   l - the region label that is desired
%
% GETOPENEPSUBSETFROMIMGMAP accepts the following parameter-value pairs
%   'repack' {true} | false
%       - Whether to repack the triangulation to remove unused vertices

% GETOPENEPSUBSETFROMIMGMAP uses labelOpenEPMapFromImgMap.m to label the
% triangulation in userdata, then creates a new triangulation containing
% only those cells which is returned in tr.
%
% TODO: Add an additional output argument which is an edited version of
% userdata containing on the data of interest; i.e. splitting userdata into
% sections.
%
% Author: Steven Williams (2021) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% 
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 3; % UPDATE VALUE
dorepack = true;
if nargin > nStandardArgs && ~isempty(varargin{1})
    for i = 1:2:nargin-1
        switch lower(varargin{i})
            case 'repack'
                dorepack = varargin{i+1};
        end
    end
end

cellLabels = labelOpenEPMapFromImgMap(userdata, mOpenEP3D);

iTriToKeep = cellLabels == l;

triangles = userdata.surface.triRep.Triangulation;
triangles(~iTriToKeep,:) = [];

X = userdata.surface.triRep.X;

tr = TriRep(triangles, X(:,1), X(:,2), X(:,3));

if dorepack
    tr = repack(tr);
end

end