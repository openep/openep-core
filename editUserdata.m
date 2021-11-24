function userdata2 = editUserdata(userdata)
% EDITUSERDATA Graphically remove regions from a Carto dataset
%
% Usage:
%   userdata2 = editUserdata(userdata)
% Where:
%   userdata  - is the original Carto dataset
%   userdata2 - is the new Carto dataset with elements removed
%
% EDITUSERDATA uses EDITTRIANGULATION to remove triangles from a TriRep
% object. Controls:
%   Left click          - select triangles to remove
%   Shift-Left click    - select triangles to keep
%   Ctrl-Left click     - select area up to the boundary
%   d                   - done
%
% Author: Steven Williams (2016) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------
%
% ---------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

tr = getMesh(userdata);
[tr2, isVertUsed] = editTriangulation(tr2);

userdata2 = userdata;
userdata2.surface.triRep = tr;

userdata2.surface.isVertexAtRim(~isVertUsed) = [];
userdata2.surface.act_bip(~isVertUsed,:) = [];

end