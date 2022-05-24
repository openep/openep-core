function userdata = import_ensitex(varargin)
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

% Create an empty OpenEP data structure
userdata = openep_createuserdata;

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

% % next load the map file and the wave files into memory
% cMapDir = [varargin{1} filesep() 'Contact_Mapping'];
% mappingFiles = nameFiles(cMapDir, 'showhiddenfiles', false, 'extension', '.csv');
% for i = 1:numel(mappingFiles)
%     [info varnames data] = loadensitex_dxldata([varargin{1} filesep() 'Contact_Mapping' filesep() mappingFiles{i}]);
%     dataFile{i}.info = info;
%     dataFile{i}.varnames = varnames;
%     dataFile{i}.data = data;
% end

load('/Users/steven/Desktop/dataFiles.mat');

% then map the data into the OpenEP data format

% I think we have to decide whether to import an 'along' map or an 'across'
% map and these might need combining retrospectively for particular use
% cases, but could be considered as separate OpenEP data structures at
% import time. The default should possibly be to import to the 'along' as 
% this is 'along the spline' and consistent with usual practice in EP, 
% whereas 'across' may be influenced by HD grid geometry.

% dataFile{1} Map_CV_omni.csv
% dataFile{2} Wave_bi_across.csv
% dataFile{3} Wave_bi_along.csv
% dataFile{4} Wave_refs.csv
% dataFile{5} Wave_rov.csv
% dataFile{6} Wave_uni_across.csv
% dataFile{7} Wave_uni_along.csv
% dataFile{8} Wave_uni_corner.csv

userdata.electric.sampleFrequency = dataFile{5}.info.sampleFreq;

userdata.electric.tags                          = dataFile{1}.data(:,26);
userdata.electric.names                         = dataFile{1}.data(:,5);
userdata.electric.electrodeNames_bip            = dataFile{3}.data(:,1);
userdata.electric.egmX                          = [str2double(dataFile{1}.data(:,7)) str2double(dataFile{1}.data(:,8)) str2double(dataFile{1}.data(:,9))];
userdata.electric.egm                           = cell2mat(dataFile{3}.data(:,6));
userdata.electric.electrodeNames_uni            = dataFile{7}.data(:,1);
userdata.electric.egmUniX(:,:,1)                = [str2double(dataFile{1}.data(:,60)) str2double(dataFile{1}.data(:,61)) str2double(dataFile{1}.data(:,62))];
userdata.electric.egmUniX(:,:,2)                = [str2double(dataFile{1}.data(:,64)) str2double(dataFile{1}.data(:,65)) str2double(dataFile{1}.data(:,66))];
userdata.electric.egmUni(:,:,1)                 = cell2mat(dataFile{8}.data(:,6));
userdata.electric.egmUni(:,:,2)                 = cell2mat(dataFile{7}.data(:,6));
% Hard coding a channel just now; but ask 'which signal is a good reference?' of the user
userdata.electric.egmRef                        = userdata.electric.egm; % only for now
% Hard coding a channel just now; but ask 'which signal is a good reference?' of the user
userdata.electric.ecg                           = userdata.electric.egm; % only for now
userdata.electric.annotations.woi               = round([str2double(dataFile{1}.data(:,23)) str2double(dataFile{1}.data(:,24))] ./ 1000 .* userdata.electric.sampleFrequency);
userdata.electric.annotations.referenceAnnot    = str2double(dataFile{1}.data(:,74));
userdata.electric.annotations.mapAnnot          = round(str2double(dataFile{3}.data(:,5)) ./ userdata.electric.sampleFrequency .* userdata.electric.sampleFrequency);
userdata.electric.voltages.bipolar              = str2double(dataFile{1}.data(:,32));
userdata.electric.votlages.unipolar             = str2double(dataFile{1}.data(:,58));
userdata.electric.egmSurfX                      = [str2double(dataFile{1}.data(:,10)) str2double(dataFile{1}.data(:,11)) str2double(dataFile{1}.data(:,12))];
userdata.electric.barDirection                  = [str2double(dataFile{1}.data(:,13)) str2double(dataFile{1}.data(:,14)) str2double(dataFile{1}.data(:,15))];
userdata.electric.include                       = str2double(dataFile{1}.data(:,17));




% infoMapping = { ...
%     'electric.sampleFrequency'                'Wave_rov.csv'                              'sampleFreq' ...
%     ; ...
% };
% 
dataMapping = { ...
       'electric.tags'                        'Map_CV_omni.csv'                           'annot' ...
    ;  'electric.names'                       'Map_CV_omni.csv'                           '(Point #)' ...
    ;  'electric.electrodeNames_bip'          'Wave_bi_along.csv'                         'Trace' ...
    ;  'electric.egmX'                        'Map_CV_omni.csv'                           'roving x, roving y, roving z' ...
    ;  'electric.egm'                         'Wave_bi_along'                             'signals' ...
    ;  'electric.electrodeNames_uni'          'Wave_uni_along'                            'Trace' ...
    ;  'electric.egmUniX'                     'Map_CV_omni.csv'                           'Uni_CornerX, Uni_CornerY, Uni_CornerZ; Uni_AlongX, Uni_AlongY, Uni_AlongZ' ...
    ;  'electric.egmUni'                      'Wave_uni_corner.csv; Wave_uni_along.csv'   'signals' ...
    ;  'electric.egmRef'                      'Wave_refs.csv'                             'signals' ... % Ask, 'which signal is a good reference'
    ;  'electric.ecg'                         'Wave_refs.csv'                             'signals' ... % Ask, 'which signal is a good ECG'
    ;  'electric.annotations.woi'             'Map_CV_omni.csv'                           'left curtain (ms), right curtain (ms)' ... % this is in ms relative to samples!
    ;  'electric.annotations.referenceAnnot'  'Map_CV_omni.csv'                           'RefTick' ... % this is in samples!
    ;  'electric.annotations.mapAnnot'        'Wave_bi_along.csv'                         'rovTime (wave samples)' ... % this is _presumably_ in samples!
    ;  'electric.voltages.bipolar'            'Map_CV_omni.csv'                           'pp_Valong' ...
    ;  'electric.voltages.unipolar'           'Map_CV_omni.csv'                           'unipoleMaxPP' ...
    ;  'electric.impedances.time'             ''                                          '' ... % We do not seem to have impedance data
    ;  'electric.impedances.value'            ''                                          '' ... % We do not seem to have impedance data
    ;  'electric.egmSurfX'                    'Map_CV_omni.csv'                           'surface x, surface y, surface z' ...
    ;  'electric.barDirection'                'Map_CV_omni.csv'                           'normal x, normal y, normal z' ...
    ;  'electric.include'                     'Map_CV_omni.csv'                           'utilized' ...
    };
% 
% % Deal with the dataMapping array
% for i = 1:size(dataMapping)
% 
%     fieldNames = strsplit(dataMapping{i,1}, '.');
% 
%     % parse the fieldnames
%     switch numel(fieldNames)
%         case 1
%             userdata.(fieldNames{1}) = dataMapping{i,1};
%         case 2
%             userdata.(fieldNames{1}).(fieldNames{2}) = dataMapping{i,1};
%         case 3
%             userdata.(fieldNames{1}).(fieldNames{2}).(fieldNames{3}) = dataMapping{i,1};
%         case 4
%             userdata.(fieldNames{1}).(fieldNames{2}).(fieldNames{3}).(fieldNames{4}) = dataMapping{i,1};
%         otherwise
%             error('OPENEP/IMPORT_ENSITEX: Code not yet implemented for more than 4 sub fields')
%     end
% 
% end
% 
% % Deal wtih the infoMapping array
% for i = 1:size(infoMapping)
% 
%     fieldNames = strsplit(infoMapping{i,1}, '.');
% 
%     % parse the fieldnames
%     switch numel(fieldNames)
%         case 1
%             userdata.(fieldNames{1}) = infoMapping{i,1};
%         case 2
%             userdata.(fieldNames{1}).(fieldNames{2}) = infoMapping{i,1};
%         case 3
%             userdata.(fieldNames{1}).(fieldNames{2}).(fieldNames{3}) = infoMapping{i,1};
%         case 4
%             userdata.(fieldNames{1}).(fieldNames{2}).(fieldNames{3}).(fieldNames{4}) = infoMapping{i,1};
%         otherwise
%             error('OPENEP/IMPORT_ENSITEX: Code not yet implemented for more than 4 sub fields')
%     end
% 
% end

 
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