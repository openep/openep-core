function data = getSurfaceData( userdata, datatype )
% GETSURFACEDATA Returns surface mapping data from userdata
%
% Usage:
%   data = getSurfaceData( userdata, datatype )
% Where:
%   userdata  - see importcarto_mem
%   datatype  - the required data. Must be one of:
%       'act', 'bip', 'uni', 'imp', 'frc'
%   data      - The returned surface mapping data
%
% GETSURFACEDATA Returns surface mapping data from userdata. Data type is
% specified by the 'datatype' argument:
%   'act' - activation time
%   'bip' - bipolar voltage
%   'uni' - unipolar voltage
%   'imp' - impedance
%   'frc' - contact force
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% data = getSurfaceData(userdata, 'bip');
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

switch datatype
    case 'act'
        data = userdata.surface.act_bip(:,1);
    case 'bip'
        data = userdata.surface.act_bip(:,2);
    case 'uni'
        data = userdata.surface.uni_imp_frc(:,1);
    case 'imp'
        data = userdata.surface.uni_imp_frc(:,2);
    case 'frc'
        data = userdata.surface.uni_imp_frc(:,3);
end
