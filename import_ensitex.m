function [userdata] = import_ensitex(varargin)
% IMPORT_ENSITEX is used to import an EnsiteX case.
%
% Usage:
%   userdata = import_ensitex()
%   userdata = import_ensitex(directory)
%   userdata = import_ensitex( ... , Name,Value ... )
%
% Where:
%   directory - an absolute folder path (if empty, user will be asked)
%   userdata - an OpenEP data structure
%
% Author: Steven Williams (2022)
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

% General data
userdata.systemName = 'ensitex';
userdata.notes{1} = [date() ': Created'];
userdata.ensiteXFolder = varargin{1};

% Load the model groups
info = loadprecision_modelgroups([varargin{1} filesep() 'Contact_Mapping_Model.xml']);

% assign the geometry data
if isfield(info.dxgeo, 'triangles') && isfield(info.dxgeo, 'vertices')
    TRI = info.dxgeo.triangles;
    X = info.dxgeo.vertices(:,1);
    Y = info.dxgeo.vertices(:,2);
    Z = info.dxgeo.vertices(:,3);
    userdata.surface.triRep = TriRep(TRI, X, Y, Z);
    FF = freeBoundary(userdata.surface.triRep);
    isVertexAtRim = false(size(userdata.surface.triRep.X,1),1);
    if ~isempty(FF)
        isVertexAtRim(FF(:,1)) = true;
    end
    userdata.surface.isVertexAtRim = isVertexAtRim;
end

% deal with labels
% for the test case this is 0, 1 or 2
allLabels = info.dxgeo.surface_of_origin;
labels = unique(allLabels);
cMap = colormap(parula(numel(labels)));
faceColors = trFaceToVertData(userdata.surface.triRep, allLabels);
hSurf = drawMap(userdata, 'type', 'none');

colorShell(hSurf, [X Y Z], faceColors, Inf ...
    , 'showcolorbar', 'show' ...
    , 'coloraxis', [0 numel(labels)] ...
    , 'interpolation', 'off' ...
    , 'usrColorMap', cMap ...
    , 'datatype', 'labels' ...
    );

% next load the map file

% then load the wave files






% 
% if exist('dxldata', 'var')
%     if isfield(dxldata, 'sampleFreq')
%         userdata.electric.sampleFrequency = dxldata.sampleFreq;
%     end
% end
% 
% % Geometry data - NOT COMPLETE
% 
% 
% % iMap = 1;
% % for i = 1:numel(data.modelgroups)
% %     for j = 1:numel(data.modelgroups(i).dxgeo)
% %         TRI = data.modelgroups(i).dxgeo(j).triangles;
% %         X = data.modelgroups(i).dxgeo(j).vertices(:,1);
% %         Y = data.modelgroups(i).dxgeo(j).vertices(:,2);
% %         Z = data.modelgroups(i).dxgeo(j).vertices(:,3);
% %         tr{iMap} = TriRep(TRI, X, Y, Z);
% %         iMap = iMap + 1;
% %     end
% % end
% 
% 
% 
% 
% % Read the data
% if nargin == 0
%     data = importprecision();
% else
%     data = importprecision('direc', varargin{1});
% end
% 
% % Load DxL data
% dxldata = importprecision_dxldata(data.directory);
% 
% % Create a blank OpenEP data structure
% userdata = openep_createuserdata();
% 
% 
% 
% 
% 
% % find the Dx Landmark Geo file index
% iDxInd = [];
% for i = 1: numel(data.modelgroups)
%     if strstartcmpi('St. Jude Medical Dx Landmark Geo', data.modelgroups(i).dxgeo.comment{i})
%         iDxInd = i;
%         break; % after finding the Dx Landmark Geo file we can proceed
%     end
% end
% if isempty(iDxInd)
%     error('OPENEP/IMPORT_PRECISION: No DxLandmarkGeo.xml file located.')
% end
% TRI = data.modelgroups(iDxInd).dxgeo.triangles;
% X = data.modelgroups(iDxInd).dxgeo.vertices(:,1);
% Y = data.modelgroups(iDxInd).dxgeo.vertices(:,2);
% Z = data.modelgroups(iDxInd).dxgeo.vertices(:,3);
% tr = TriRep(TRI, X, Y, Z);
% userdata = setMesh(userdata, tr);
% surfaceData = data.modelgroups(iDxInd).dxgeo.surface_of_origin;
% userdata = setSurfaceProperty(userdata, 'name', 'surfaceOfOrigin', 'map', surfaceData, 'definedOn', 'elements');
% 
% 
% % Surface data - TODO
% % THERE DOES NOT APPEAR TO BE VOLTAGE MAPPING DATA AVAILABLE BUT THERE DOES
% % APPEAR TO BE ACTIVATION TIME DATA AVAILABLE IN THE Dx Landmark Geo data file
% % Get the activation time 
% if isfield(data.modelgroups(iDxInd).dxgeo, 'act')
%     act = data.modelgroups(iDxInd).dxgeo.act;
% else 
%     act = repmat(NaN, size(X));
% end
% if isfield(data.modelgroups(iDxInd).dxgeo, 'bip')
%     bip = data.modelgroups(iDxInd).dxgeo.bip;
% else
%     bip = repmat(NaN, size(X));
% end
% userdata.surface.act_bip = [act bip];
% % userdata.surface.uni_imp_frc = 
% 
% % if isfield(data.modelgroups(i).dxgeo(j).vertices(:,3);)
% % userdata = setSurfaceProperty(userdata, 'name', 'surfaceregion', info.dxgeo.surface_of_origin);
% 
% 
% 
% 
% % TODO: HOW TO DEAL WITH GEOMETRY *** Github Issue #42: https://github.com/openep/openep-core/issues/42 ***
% % (1) Find the dxgeo with the first comment 'St. Jude Medical Dx Landmark Geo data export; file format revision 0'
% %     - best to do this by identifying the comment which contains the string 'Dx Landmark Geo'
% % (2) Use the geometry described in this file as userdata.surface.triRep
% % (3) Find the dxgeo with the comment 'St. Jude Medical Model Groups data export; file format revision 0'
% %     - best to do this be identifying the comment which contains the string 'Model Groups'
% % (4) Use this file to tag each polygon in userdata.surface.triRep as
% % belongining to one or more geometries. Will need to think of how/where to
% % store this information, but it probably needs another data field and may
% % be specific to Precision
% % (5) Note that some dxgeo's do not contain a comment. These *probably*
% % refer to image data sets merged into the system. Will need to think of
% % how/where to store this information, but it probably needs another data
% % field, and is not likley to be specific to Precision
% 
% % this section deals with geometry
% % sometimes there may be model groups; sometimes not
% 
% 
% 
% 
% 
% 
% % Electric data - PARTIALLY COMPLETE
% % userdata.electric.tags = ;
% % userdata.electric.names = ;
% if length(dxldata) > 2
%     warning(['Currently unable to process more than one combined set ',...
%         'of experiments (one unipolar and one bipolar'])
% end
% for i_dxl = 1:length(dxldata)
%     if dxldata(1).bipole
%         userdata.electric.electrodeNames_bip = dxldata(i_dxl).rovtrace_pts';
%         userdata.electric.egmX = [dxldata(i_dxl).rovingx',...
%             dxldata(i_dxl).rovingy', dxldata(i_dxl).rovingz'];
%         userdata.electric.egmSurfX = [dxldata(i_dxl).surfPtx',...
%             dxldata(i_dxl).surfPty' dxldata(i_dxl).surfPtz'];
%         userdata.electric.egmRef = dxldata(i_dxl).rovtrace'; % TODO: import the reference egm
%         userdata.electric.egm = dxldata(i_dxl).rovtrace';
%         userdata.electric.annotations.referenceAnnot = dxldata(i_dxl).refLAT';
%         userdata.electric.annotations.mapAnnot = dxldata(i_dxl).rovLAT';
%         userdata.electric.annotations.woi = -userdata.electric.annotations.referenceAnnot;
%         userdata.electric.annotations.woi(:,2) = length(userdata.electric.egm)-userdata.electric.annotations.referenceAnnot;
%         userdata.electric.voltages.bipolar = dxldata(i_dxl).peak2peak';
%     else
%         userdata.electric.electrodeNames_uni = dxldata(i_dxl).rovtrace_pts';
%         userdata.electric.egmUniX = [dxldata(i_dxl).rovingx',...
%             dxldata(i_dxl).rovingy', dxldata(i_dxl).rovingz'];
%         userdata.electric.egmUniSurfX = [dxldata(i_dxl).surfPtx',...
%             dxldata(i_dxl).surfPty' dxldata(i_dxl).surfPtz'];
%         userdata.electric.egmUniRef = dxldata(i_dxl).rovtrace';
%         userdata.electric.egmUni = dxldata(i_dxl).rovtrace';
%         userdata.electric.annotations.referenceAnnotUni = dxldata(i_dxl).refLAT';
%         userdata.electric.annotations.mapAnnotUni = dxldata(i_dxl).rovLAT';
%         userdata.electric.voltages.unipolar = dxldata(i_dxl).peak2peak';
%     end
% end
% % userdata.electric.egm = rovtrace';
% % TODO: IMPORTING UNIPOLE DATA *** Github Issue #43: https://github.com/openep/openep-core/issues/43***
% % userdata.electric.electrodeNames_uni = ; 
% % userdata.electric.egmUniX = ;
% % userdata.electric.egmUni = ;
% % userdata.electric.ecg =
% % userdata.electric.annotations.woi = 
% % userdata.electric.voltages.unipolar = 
% % userdata.electric.impedances.time = 
% % userdata.electric.impedances.value = 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % Ablation data - TODO
% % userdata.rf.originaldata.force.time = 
% % userdata.rf.originaldata.force.force =
% % userdata.rf.originaldata.force.axialangle = 
% % userdata.rf.originaldata.force.lateralangle = 
% % userdata.rf.originaldata.force.position = 
% % userdata.rf.originaldata.ablparams.time = 
% % userdata.rf.originaldata.ablparams.power =
% % userdata.rf.originaldata.ablparams.impedance = 
% % userdata.rf.originaldata.ablparams.distaltemp = 
% 
% end