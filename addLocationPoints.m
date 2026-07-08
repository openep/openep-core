function userdata = addLocationPoints(userdata, X, tags, names)
% ADDLOCATIONPOINTS Manually add location points to data
%
% Usage:
%   userdata = addLocationPoints(userdata, X, tags, names)
% Where:
%   userdata  - the OpenEP data structure
%   X  - locations
%   tags  - tags
%   names - names
%
% addLocationPoints does not accept parameter value pairs
%
% ADDLOCATIONPOINTS Manually add location points to data
%
% Author: Steven Williams (2024)
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% % copy and paste co-ordinates from a CSV file
% tags = cell(size(X,1), 1);
% tags(:) = {'ablation'};
% names = cell(size(X,1), 1);
% names(:) = {'abl'};
% userdata = addLocationPoints(userdata, X, tags, names)
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

numPointsToAdd = size(X, 1);

for i = 1:numPointsToAdd
userdata.electric.tags{end+1} = tags{i};
userdata.electric.names{end+1} = [names{i} num2str(i)];
userdata.electric.egmX(end+1,:) = X(i,:);
userdata.electric.egm(end+1,:) = NaN;
userdata.electric.egmRef(end+1,:) = NaN(1,1000);
userdata.electric.annotations.woi(end+1,:) = [NaN NaN];
userdata.electric.annotations.referenceAnnot(end+1,:) = NaN;
userdata.electric.annotations.mapAnnot(end+1,:) = NaN;
userdata.electric.voltages.bipolar(end+1,:) = NaN;

tr = getMesh(userdata, 'type', 'triangulation');
surfX = findclosestvertex(tr, X(i,:), true);
userdata.electric.egmSurfX(end+1,:) = userdata.surface.triRep.X(surfX,:);
userdata.electric.LATs(end+1,:) = NaN;
userdata.electric.electrodeNames_uni(end+1,:) = NaN;
userdata.electric.include(end+1) = false;
userdata.electric.discarded(end+1) = true;


end