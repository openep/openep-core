function [ userdata_ds ] = downSampleDataset( userdata, sampleDensity )
% DOWNSAMPLEDATASET Returns a new OpenEP dataset at the required sampling
% density
%
% Usage:
%   [ userdata_ds ] = DOWNSAMPLEDATASET( userdata, sampleDensity )
% Where:
%   userdata  - the input
%   userdata_ds  - the output
%
% DOWNSAMPLEDATASET Does not check the distribution of points but just
% returns the required number of points per map
%
% Author: Steven Williams (2021)
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

% Work out the surface area and original density
surfaceArea = getArea(userdata);
numPts = getNumPts(userdata);
densityOriginal = numPts/surfaceArea;
disp(['Original density = ' num2str(densityOriginal) 'points/cm2'])

% Calculate and identify the required number and points
requiredNumberOfPoints = ceil(surfaceArea * sampleDensity);
if requiredNumberOfPoints > numPts
    error('OPENEP/DOWNSAMPLEDATASET: Original sampling density is insufficient to fulfill this request');
end
indices = randperm(numPts);
iKeep = indices(1:requiredNumberOfPoints);
iRemove = indices(requiredNumberOfPoints+1:end);

% Remove the data we do not need
userdata_ds = userdata;
userdata_ds.electric.tags(iRemove,:) = [];
userdata_ds.electric.names(iRemove,:) = [];
userdata_ds.electric.electrodeNames_bip(iRemove,:) = [];
userdata_ds.electric.egmX(iRemove,:) = [];
userdata_ds.electric.egm(iRemove,:) = [];
userdata_ds.electric.electrodeNames_uni(iRemove,:) = [];
userdata_ds.electric.egmUniX(iRemove,:,:) = [];
userdata_ds.electric.egmUni(iRemove,:,:) = [];
userdata_ds.electric.egmRef(iRemove,:) = [];
userdata_ds.electric.ecg(iRemove,:) = [];
userdata_ds.electric.annotations.woi(iRemove,:) = [];
userdata_ds.electric.annotations.referenceAnnot(iRemove,:) = [];
userdata_ds.electric.annotations.mapAnnot(iRemove,:) = [];
userdata_ds.electric.voltages.bipolar(iRemove,:) = [];
userdata_ds.electric.voltages.unipolar(iRemove,:) = [];
userdata_ds.electric.impedances.time(:,iRemove) = [];
userdata_ds.electric.impedances.value(:,iRemove) = [];
userdata_ds.electric.egmSurfX(iRemove,:) = [];
userdata_ds.electric.barDirection(iRemove,:) = [];

% Work out the new density
numPts_ds = getNumPts(userdata_ds);
density_ds = numPts_ds/surfaceArea;
disp(['New density = ' num2str(density_ds) 'points/cm2'])

% Add a note
userdata_ds.notes{end+1} = [date ': downsampled']

end