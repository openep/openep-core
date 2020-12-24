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

faces = userdata.surface.triRep.Triangulation;

end
