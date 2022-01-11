function faces = getFaces( userdata )
% GETFACES Returns the faces referenced by userdata
%
% Usage:
%   faces = getFaces( userdata )
% Where:
%   userdata  - see importcarto_mem
%   faces - all the faces
%
% GETFACES Returns the faces referenced by userdata
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% faces = getFaces( userdata )
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

if isa(userdata.surface.triRep, 'TriRep')
    faces = userdata.surface.triRep.Triangulation;
elseif isa(userdata.surface.triRep, 'triangulation')
    faces = userdata.surface.triRep.ConnectivityList;
elseif isa(userdata.surface.triRep, 'struct')
    faces = userdata.surface.triRep.Triangulation;
else
    error('OPENEP/getFaces: userdata.surface.TriRep should be a TriRep, struct or a triangulation object')
end

end
