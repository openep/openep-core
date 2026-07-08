function info = loadprecision_modelgroups(filename)
% LOADPRECISION_MODELGROUPS loads a Precision file.
% Usage:
%   info = loadprecision_modelgroups(filename)
% Where:
%   filename is a full path to a file.
%   info.dxgeo is a structure containing the geometry and the labels
%   info.dataElement is 'modelgroups' - used to identify the data in
%                           compatibility with loadprecision_wavefile
%
% LOADPRECISION_MODELGROUPS ..........
% If filename does not contain a file with a valid format, info is left
% empty and a warning generated with msgID = 'LoadPrecision:InvalidFile'.

% Author: Nick Linton (2009)
% Modifications
%   Steven Williams 2017    Convert to Precision
%   Nick Linton 2017        Updated to info.dataElement and info.dxgeo in
%                           order to keep compatibility with other
%                           loadprecision_ code.
%
% Info on Code Testing:
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

info = [];

[~,~,ext] = fileparts(filename);
if ~strcmpi(ext, '.xml')
    warning('LoadPrecision:InvalidFile','LOADPRECISION_MODELGROUPS: .xml file expected');
    return
end

[tree , treeNames] = xml_read(filename);

% check version
if ~strcmpi(tree.DIFHeader.Version, 'SJM_DIF_5.0')
    warning('LoadPrecision:InvalidFile','LOADPRECISION_MODELGROUPS: unexpected version number.')
    return
end

for i = 1:tree.DIFBody.Volumes.ATTRIBUTE.number
    % pre-populate .act, .bip and .uni
    dxgeo(i).act = [];
    dxgeo(i).bip = [];
    dxgeo(i).uni = [];
    if isfield(tree.DIFBody.Volumes.Volume(i), 'Vertices')
        dxgeo(i).vertices = str2num(tree.DIFBody.Volumes.Volume(i).Vertices.CONTENT);
    end
    if isfield(tree.DIFBody.Volumes.Volume(i), 'Polygons')
        dxgeo(i).triangles = str2num(tree.DIFBody.Volumes.Volume(i).Polygons.CONTENT);
    end
    if isfield(tree.DIFBody.Volumes.Volume(i), 'Normals')
        dxgeo(i).normals = str2num(tree.DIFBody.Volumes.Volume(i).Normals.CONTENT);
    end
    if isfield(tree.DIFBody.Volumes.Volume(i), 'COMMENT')
        dxgeo(i).comment = tree.DIFBody.Volumes.Volume(i).COMMENT;
    end
    if isfield(tree.DIFBody.Volumes.Volume(i).ATTRIBUTE, 'name')
        dxgeo(i).name = tree.DIFBody.Volumes.Volume(i).ATTRIBUTE.name;
    end
    if isfield(tree.DIFBody.Volumes.Volume(i), 'Surface_of_origin')
        dxgeo(i).surface_of_origin = str2num(tree.DIFBody.Volumes.Volume(i).Surface_of_origin.CONTENT);
    end
    if isfield(tree.DIFBody.Volumes.Volume(i), 'Map_data')
        % Identify the type of data that we have. LAT calculated from a
        % bipolar map will be stored in .act; voltage calcualted from a
        % bipolar map will be stored in .bip; voltage calculated from a
        % unipolar map will be stored in .uni. Anything else will be stored
        % in .mapData with a description given in .mapType. In a future
        % version we may want to modify this logic so that all surface maps
        % are stored as mapData/mapType, to align with the concept of
        % surface properties in the OpenEP data format.
        stub = 'Data values at each vertex of DxL map';
        iComment = [];
        for j=1:numel(dxgeo.comment)
            if isstr(dxgeo.comment{j})
                if strstartcmpi(stub, dxgeo.comment{j})
                    iComment = j;
                end
            end
        end
        if isempty(iComment)
            warning(['OPENEP/LOADPRECISION_MODELGROUPS: Map data was ...' ...
                'found but no comment describing its type was identified. ...' ...
                'The data will be stored in mapdata, and the description will be set to UNKNOWN']);
            dxgeo(i).mapdata = str2num(tree.DIFBody.Volumes.Volume(i).Map_data.CONTENT);
            dxgeo(i).maptype = 'UNKNWOWN';
        else
            dataTypeString = dxgeo.comment{iComment};
            if contains(dataTypeString, 'P-P Voltage') && contains(dataTypeString, 'bi')
                dxgeo(i).bip = str2num(tree.DIFBody.Volumes.Volume(i).Map_data.CONTENT);
            elseif contains(dataTypeString, 'P-P Voltage') && contains(dataTypeString, 'uni')
                dxgeo(i).uni = str2num(tree.DIFBody.Volumes.Volume(i).Map_data.CONTENT);
            elseif contains(dataTypeString, 'LAT Isochronal') && contains(dataTypeString, 'bi')
                dxgeo(i).act = str2num(tree.DIFBody.Volumes.Volume(i).Map_data.CONTENT);
            elseif contains(dataTypeString, 'LAT Isochronal')
                warning(['OPENEP/LOADPRECISION_MODELGROUPS: Map data was  ...' ...
                    'found but its type was not identifiable as act, ...' ...
                    'bip or uni. The data will be stored in mapdata, and the description in maptype']);
                dxgeo(i).mapdata = str2num(tree.DIFBody.Volumes.Volume(i).Map_data.CONTENT);
                maptype = regexp(dataTypeString, 'map\s+(.*)$', 'tokens', 'once');
                dxgeo(i).maptype = maptype{1};
            end
        end
    end
    if isfield(tree.DIFBody.Volumes.Volume(i), 'Map_status')
        dxgeo(i).map_status = str2num(tree.DIFBody.Volumes.Volume(i).Map_status.CONTENT);
    end

    %get the labels
    if isfield(tree.DIFBody.Labels, 'Label')
        for j = 1:tree.DIFBody.Labels.ATTRIBUTE.number
            dxgeo(i).labels(j).xyz = tree.DIFBody.Labels.Label(j).CONTENT;
            dxgeo(i).labels(j).name = tree.DIFBody.Labels.Label(j).ATTRIBUTE;
        end
    end
end

info.dataElement = 'modelgroups';
info.dxgeo = dxgeo;
info.fileLoaded = filename;
end

