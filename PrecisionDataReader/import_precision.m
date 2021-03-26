function userdata = import_precision(varargin)
% IMPORT_PRECISION is used to import Precision data.
%
% Usage:
%   userdata = import_precision()
%   userdata = import_precision(directory)
%   userdata = import_precision( ... , Name,Value ... )
%
% Where:
%   directory - an absolute folder path (if empty, user will be asked)
%   userdata - see importcarto_mem
%
% IMPORT_PRECISION accepts the following parameter-value pairs
%   'filematch'     {} | String
%       A string that gives a partial match to the file(s) to be loaded.
%       e.g. {'bipol_RAW', 'Location'}. Not case sensitive and does not
%       need to be a full match. This is useful if you do not want to read
%       in all files (save time + memory).
%
% IMPORT_PRECISION is a wrapper function for importprecison.m and
% translates the output of importprecision.m into OpenEP format userdata.
%
% Author: Steven Williams / Nick Linton (2021)
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

% Read the data
data = import_precision(varargin);

% Translate into OpenEP format
userdata = openep_createuserdata();

% General data - COMPLETE
userdata.systemName = 'precision';
userdata.notes{1} = [date() ': Created'];
if isfield(data.directory)
    userdata.precisionFolder = data.directory;
end
if isfield(data.epcath_bip_raw.sampleFreq)
    userdata.electric.sampleFrequency = data.epcath_bip_raw.sampleFreq;
end

% Electric data - TODO
% userdata.electric.tags = 
% userdata.electric.names = 
% userdata.electric.electrodeNames_bip = 
% userdata.electric.egmX =
% userdata.electric.egmSurfX =
% userdata.electric.egm =
% userdata.electric.electrodeNames_uni = 
% userdata.electric.egmUniX = 
% userdata.electric.egmUni =
% userdata.electric.egmRef = 
% userdata.electric.ecg =
% 
% userdata.electric.annotations.woi = 
% userdata.electric.annotations.referenceAnnot = 
% userdata.electric.annotations.mapAnnot =
% 
% userdata.electric.voltages.bipolar = 
% userdata.electric.voltages.unipolar = 
% 
% userdata.electric.impedances.time = 
% userdata.electric.impedances.value = 

% Geometry data - COMPLETE
if isfield(data.modelgroups.dxgeo.triangles) && isfield(data.modelgroups.dxgeo.vertices)
    TRI = data.modelgroups.dxgeo.triangles;
    X = data.modelgroups.dxgeo.vertices(:,1);
    Y = data.modelgroups.dxgeo.vertices(:,2);
    Z = data.modelgroups.dxgeo.vertices(:,3);
    userdata.surface.triRep = TriRep(TRI, X, Y, Z);
    FF = freeBoundary(userdata.surface.triRep);
    isVertexAtRim = false(size(userdata.surface.triRep.X,1),1);
    isVertexAtRim(FF(:,1)) = true;
    userdata.surface.isVertexAtRim = isVertexAtRim;
end

% Surface data TODO
% userdata.surface.act_bip = 
% userdata.surface.uni_imp_frc = 

% Ablation data - TODO
% userdata.rf.originaldata.force.time = 
% userdata.rf.originaldata.force.force =
% userdata.rf.originaldata.force.axialangle = 
% userdata.rf.originaldata.force.lateralangle = 
% userdata.rf.originaldata.force.position = 
% userdata.rf.originaldata.ablparams.time = 
% userdata.rf.originaldata.ablparams.power =
% userdata.rf.originaldata.ablparams.impedance = 
% userdata.rf.originaldata.ablparams.distaltemp = 

end