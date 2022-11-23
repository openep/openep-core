function userdata = addLandMarkPoint(userdata, pointName, position)
% ADDLANDMARKPOINT Adds a landmark point named pointName with the
% co-ordinates given in position
%
% Usage:
%   userdata = addLandMarkPoint(userdata, pointName, position)
% Where:
%   userdata  - an OpenEP dataset
%   pointName  - the required point name
%   position  - the point co-ordinates
%
% ADDLANDMARKPOINT does not accept any parameter-value pairs
%
% ADDLANDMARKPOINT adds the pointName to userdata.electric.names, and the
% type 'Tissue' to userdata.electric.tags. It then adds NaNs to extend
% arrays as needed. Note this has not been implemented for datasets with
% unipolar electrogram data.
%
% Author: Steven Williams (2022)
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

userdata.electric.names = [userdata.electric.names; pointName];
userdata.electric.tags = [userdata.electric.tags; 'Tissue'];
userdata.electric.egmX = [userdata.electric.egmX; position];
userdata.electric.egm = [userdata.electric.egm; NaN(1, size(userdata.electric.egm,2))];
userdata.electric.egmRef = [userdata.electric.egmRef; NaN(1, size(userdata.electric.egmRef,2))];
userdata.electric.annotations.woi = [userdata.electric.annotations.woi; [NaN NaN]];
userdata.electric.annotations.referenceAnnot = [userdata.electric.annotations.referenceAnnot; 0];
userdata.electric.annotations.mapAnnot = [userdata.electric.annotations.mapAnnot; 0];
userdata.electric.voltages.bipolar = [userdata.electric.voltages.bipolar; NaN];
[vertices, ~] = findclosestvertex(userdata.surface.triRep.X, position);
coord = userdata.surface.triRep.X(vertices,:);
userdata.electric.egmSurfX = [userdata.electric.egmSurfX; coord];
userdata.electric.electrodeNames_uni = [userdata.electric.electrodeNames_uni; [NaN NaN]];
userdata.electric.include = [userdata.electric.include; 0];

