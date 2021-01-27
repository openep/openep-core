function iPoint = getMappingPointsWithinWoI( userdata )
% GETMAPPINGPOINTSWITHINWOI Returns the indices of the mapping points with
% annotated local activation time within the window of interest
%
% Usage:
%   iPoint = getMappingPointsWithinWoI( userdata )
% Where:
%   userdata  - see importcarto_mem
%   iPoint  - logicla array list of valid points; indexes into userdata.electric
%
% GETMAPPINGPOINTSWITHINWOI Returns the indices of the mapping points with
% annotated local activation time within the window of interest
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% See also GETMAPPINGPOINTSWITHINMESH
%
% Info on Code Testing:
% ---------------------------------------------------------------
% iPoint = getMappingPointsWithinWoI( userdata );
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

referenceAnnot = userdata.electric.annotations.referenceAnnot;
woi = userdata.electric.annotations.woi;
woi = woi + [referenceAnnot referenceAnnot];
mapAnnot = userdata.electric.annotations.mapAnnot;
iPoint = true(size(mapAnnot));
for i = 1:numel(mapAnnot)
    if mapAnnot(i)<woi(i,1) || mapAnnot(i)>woi(i,2)
        iPoint(i) = false;
    end
end

end