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
if isfield(data.epcath_bip_raw, 'sampleFreq')
    userdata.electric.sampleFrequency = data.epcath_bip_raw.sampleFreq;
end

% Geometry data - NOT COMPLETE
% if isfield(data.modelgroups.dxgeo, 'triangles') && isfield(data.modelgroups.dxgeo, 'vertices')
%     TRI = data.modelgroups.dxgeo.triangles;
%     X = data.modelgroups.dxgeo.vertices(:,1);
%     Y = data.modelgroups.dxgeo.vertices(:,2);
%     Z = data.modelgroups.dxgeo.vertices(:,3);
%     userdata.surface.triRep = TriRep(TRI, X, Y, Z);
%     FF = freeBoundary(userdata.surface.triRep);
%     isVertexAtRim = false(size(userdata.surface.triRep.X,1),1);
%     isVertexAtRim(FF(:,1)) = true;
%     userdata.surface.isVertexAtRim = isVertexAtRim;
% end
iMap = 1
for i = 1:numel(data.modelgroups)
    for j = 1:numel(data.modelgroups(i).dxgeo)
        TRI = data.modelgroups(i).dxgeo(j).triangles;
        X = data.modelgroups(i).dxgeo(j).vertices(:,1);
        Y = data.modelgroups(i).dxgeo(j).vertices(:,2);
        Z = data.modelgroups(i).dxgeo(j).vertices(:,3);
        tr{iMap} = TriRep(TRI, X, Y, Z);
        iMap = iMap + 1;
    end
end
% TODO HOW TO DEAL WITH GEOMETRY
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
% userdata.electric.tags = ;
% userdata.electric.names = ;
userdata.electric.electrodeNames_bip = dxldata.rovtrace_pts';
userdata.electric.egmX = [dxldata.rovingx' dxldata.rovingy' dxldata.rovingz'];
userdata.electric.egmSurfX = [dxldata.surfPtx' dxldata.surfPty' dxldata.surfPtz'];
userdata.electric.egm = rovtrace';
% userdata.electric.electrodeNames_uni = ; 
% userdata.electric.egmUniX = ;
% userdata.electric.egmUni = ;
userdata.electric.egmRef = dxldata.rovtrace';
% userdata.electric.ecg =
% userdata.electric.annotations.woi = 
userdata.electric.annotations.referenceAnnot = dxldata.refLAT';
userdata.electric.annotations.mapAnnot = dxldata.rovLAT';
userdata.electric.voltages.bipolar = dxldata.peak2peak';
% userdata.electric.voltages.unipolar = 
% userdata.electric.impedances.time = 
% userdata.electric.impedances.value = 

% Surface data - TODO
% NOTE THE SURFACE DATA WILL COME FROM THE DXL FILES I THINK
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