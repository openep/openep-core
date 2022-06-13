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
if nargin == 0
    data = importprecision();
else
    data = importprecision('direc', varargin{1});
end

% Load DxL data
dxldata = importprecision_dxldata(data.directory);

% Create a blank OpenEP data structure
userdata = openep_createuserdata();

% General data - COMPLETE
userdata.systemName = 'precision';
userdata.notes{1} = [date() ': Created'];
if isfield(data, 'directory')
    userdata.precisionFolder = data.directory;
end
if exist('dxldata', 'var')
    if isfield(dxldata, 'sampleFreq')
        userdata.electric.sampleFrequency = dxldata.sampleFreq;
    end
end

% find the Dx Landmark Geo file index
iDxInd = [];
for i = 1: numel(data.modelgroups)
    if strstartcmpi('St. Jude Medical Dx Landmark Geo', data.modelgroups(i).dxgeo.comment{i})
        iDxInd = i;
        break; % after finding the Dx Landmark Geo file we can proceed
    end
end
if isempty(iDxInd)
    error('OPENEP/IMPORT_PRECISION: No DxLandmarkGeo.xml file located.')
end
TRI = data.modelgroups(iDxInd).dxgeo.triangles;
X = data.modelgroups(iDxInd).dxgeo.vertices(:,1);
Y = data.modelgroups(iDxInd).dxgeo.vertices(:,2);
Z = data.modelgroups(iDxInd).dxgeo.vertices(:,3);
tr = TriRep(TRI, X, Y, Z);

t.X = tr.X;
t.Triangulation = tr.Triangulation;
userdata.surface.triRep = t;

surfaceData = data.modelgroups(iDxInd).dxgeo.surface_of_origin;
userdata = setSurfaceProperty(userdata, 'name', 'surfaceOfOrigin', 'map', surfaceData, 'definedOn', 'elements');


% Surface data - PARTIALLY COMPLETE
% There appears to be only activation OR voltage data in the mapping file and not both
if isfield(data.modelgroups(iDxInd).dxgeo, 'act')
    act = data.modelgroups(iDxInd).dxgeo.act;
else 
    act = repmat(NaN, size(X));
end
if isfield(data.modelgroups(iDxInd).dxgeo, 'bip')
    bip = data.modelgroups(iDxInd).dxgeo.bip;
else
    bip = repmat(NaN, size(X));
end
userdata.surface.act_bip = [act bip];
% userdata.surface.uni_imp_frc = 

% if isfield(data.modelgroups(i).dxgeo(j).vertices(:,3);)
% userdata = setSurfaceProperty(userdata, 'name', 'surfaceregion', info.dxgeo.surface_of_origin);

% TODO: HOW TO DEAL WITH GEOMETRY *** Github Issue #42: https://github.com/openep/openep-core/issues/42 ***
% (1) Find the dxgeo with the first comment 'St. Jude Medical Dx Landmark Geo data export; file format revision 0'
%     - best to do this by identifying the comment which contains the string 'Dx Landmark Geo'
% (2) Use the geometry described in this file as userdata.surface.triRep
% (3) Find the dxgeo with the comment 'St. Jude Medical Model Groups data export; file format revision 0'
%     - best to do this be identifying the comment which contains the string 'Model Groups'
% (4) Use this file to tag each polygon in userdata.surface.triRep as
% belongining to one or more geometries. Will need to think of how/where to
% store this information, but it probably needs another data field and may
% be specific to Precision
% (5) Note that some dxgeo's do not contain a comment. These *probably*
% refer to image data sets merged into the system. Will need to think of
% how/where to store this information, but it probably needs another data
% field, and is not likley to be specific to Precision

% Electric data - PARTIALLY COMPLETE
if length(dxldata) > 2
    warning(['Currently unable to process more than one combined set ',...
        'of experiments (one unipolar and one bipolar'])
end
for i_dxl = 1:length(dxldata)
    if dxldata(1).bipole
        userdata.electric.electrodeNames_bip = dxldata(i_dxl).rovtrace_pts';
        userdata.electric.egmX = [dxldata(i_dxl).rovingx',...
            dxldata(i_dxl).rovingy', dxldata(i_dxl).rovingz'];
        userdata.electric.egmSurfX = [dxldata(i_dxl).surfPtx',...
            dxldata(i_dxl).surfPty' dxldata(i_dxl).surfPtz'];
        userdata.electric.egmRef = dxldata(i_dxl).rovtrace'; % TODO: import the reference egm
        userdata.electric.egm = dxldata(i_dxl).rovtrace';
        userdata.electric.annotations.referenceAnnot = dxldata(i_dxl).refLAT'; % this is in samples
        userdata.electric.annotations.mapAnnot = dxldata(i_dxl).rovLAT'; % this is in samples

        userdata.electric.annotations.woi = 1 - userdata.electric.annotations.referenceAnnot;
        userdata.electric.annotations.woi(:,2) = size(userdata.electric.egm,2) - userdata.electric.annotations.referenceAnnot;

        userdata.electric.voltages.bipolar = dxldata(i_dxl).peak2peak';
    else
        warning('OPENEP/IMPORT_PRECISION: This code is not fully tested and likely to yield errors')
        userdata.electric.electrodeNames_uni = dxldata(i_dxl).rovtrace_pts';
        userdata.electric.egmUniX = [dxldata(i_dxl).rovingx',...
            dxldata(i_dxl).rovingy', dxldata(i_dxl).rovingz'];
        userdata.electric.egmUniSurfX = [dxldata(i_dxl).surfPtx',...
            dxldata(i_dxl).surfPty' dxldata(i_dxl).surfPtz'];
        userdata.electric.egmUniRef = dxldata(i_dxl).rovtrace';
        userdata.electric.egmUni = dxldata(i_dxl).rovtrace';
        userdata.electric.annotations.referenceAnnotUni = dxldata(i_dxl).refLAT';
        userdata.electric.annotations.mapAnnotUni = dxldata(i_dxl).rovLAT';
        userdata.electric.voltages.unipolar = dxldata(i_dxl).peak2peak';
    end
end

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