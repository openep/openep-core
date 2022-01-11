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
if ~isempty(userdata.electric.tags)
    userdata_ds.electric.tags(iRemove,:) = [];
end
if ~isempty(userdata.electric.names)
    userdata_ds.electric.names(iRemove,:) = [];
end
if ~isempty(userdata_ds.electric.electrodeNames_bip)
    userdata_ds.electric.electrodeNames_bip(iRemove,:) = [];
end
if ~isempty(userdata_ds.electric.egmX)
    userdata_ds.electric.egmX(iRemove,:) = [];
end
if ~isempty(userdata_ds.electric.egm)
    userdata_ds.electric.egm(iRemove,:) = [];
end
if ~isempty(userdata_ds.electric.electrodeNames_uni)
    userdata_ds.electric.electrodeNames_uni(iRemove,:) = [];
end
if ~isempty(userdata_ds.electric.egmUniX)
    userdata_ds.electric.egmUniX(iRemove,:,:) = [];
end
if ~isempty(userdata_ds.electric.egmUni)
    userdata_ds.electric.egmUni(iRemove,:,:) = [];
end
if ~isempty(userdata_ds.electric.egmRef)
    userdata_ds.electric.egmRef(iRemove,:) = [];
end
if ~isempty(userdata_ds.electric.ecg)
    userdata_ds.electric.ecg(iRemove,:) = [];
end
if ~isempty(userdata_ds.electric.annotations.woi)
    userdata_ds.electric.annotations.woi(iRemove,:) = [];
end
if ~isempty(userdata_ds.electric.annotations.referenceAnnot)
    userdata_ds.electric.annotations.referenceAnnot(iRemove,:) = [];
end
if ~isempty(userdata_ds.electric.annotations.mapAnnot)
    userdata_ds.electric.annotations.mapAnnot(iRemove,:) = [];
end
if ~isempty(userdata_ds.electric.voltages.bipolar)
    userdata_ds.electric.voltages.bipolar(iRemove,:) = [];
end
if ~isempty(userdata_ds.electric.voltages.unipolar)
    userdata_ds.electric.voltages.unipolar(iRemove,:) = [];
end
if ~isempty(userdata_ds.electric.impedances.time)
    userdata_ds.electric.impedances.time(:,iRemove) = [];
end
if ~isempty(userdata_ds.electric.impedances.value)
    userdata_ds.electric.impedances.value(:,iRemove) = [];
end
if ~isempty(userdata_ds.electric.egmSurfX)
    userdata_ds.electric.egmSurfX(iRemove,:) = [];
end
if ~isempty(userdata_ds.electric.barDirection)
    userdata_ds.electric.barDirection(iRemove,:) = [];
end

% Work out the new density
numPts_ds = getNumPts(userdata_ds);
density_ds = numPts_ds/surfaceArea;
disp(['New density = ' num2str(density_ds) 'points/cm2'])

% Add a note
userdata_ds.notes{end+1} = [date ': downsampled']

end