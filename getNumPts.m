function numpts = getNumPts( userdata )
% GETNUMPTS Returns the number of mapping points in the OpenEP dataset
%
% Usage:
%   numpts = getNumPts( userdata )
% Where:
%   userdata  - see importcarto_mem
%   numpts  - the number of mapping points
%
% GETNUMPTS Returns the number of mapping points in the OpenEP datasest
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
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

numpts = numel(userdata.electric.names);

end