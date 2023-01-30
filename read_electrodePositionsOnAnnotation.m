function [electrodePositions, namesRead ] = read_electrodePositionsOnAnnotation(names, pointFileName)
% READ_ELECTRODEPOSITIONSONANNOTATION gets positions of named electrodes.
% Usage:
%   [ electrogramname_bip, electrogramname_uni ] = getpointelectrogramname( point_xyz, pointFileName )
% Where:
%   names   - name of a...
%                 single unipole channel eg. '20A_1'
%                 bipole channel, eg. '20A_1-2'
%           - a cell array combining any of the above
%           - name of a SINGLE connector - e.g any of
%                          {'CS_CONNECTOR'; ...
%                           'MAGNETIC_20_POLE_A_CONNECTOR'; ...
%                           'MAGNETIC_20_POLE_B_CONNECTOR'; ...
%                           'NAVISTAR_CONNECTOR'; ...
%                           'MEC'; ...
%                           'MCC_DX_CONNECTOR'};
%           - note that if ECG names are included then NAN will be returned
%             for these
%   pointFileName - must be the full path and it is assumed that the other
%                   relevant files are in the same directory. e.g.
%                   E:\Export_PAF-01_20_2023-14-13-23\1-1-LA_P1_Point_Export.xml
%
%   electrodePositions - a double array of positions, size [numel(namesRead),3]
%   namesRead - usually the same as names (NaN is returned in the case that
%   no match is found for allowed
%
% READ_ELECTRODEPOSITIONSONANNOTATION looks in the same directory as the
% directory in which POINTFILENAME is placed. The electrogram NAMES are
% interrogated to work out the relevant CONNECTOR and that file is opened
% and read.
% Note that the bipole position is calculated as the average of the two
% electrode positions (this is not standard in Carto).
% If the names of ECG signals are passed then NaN values is returned.
% If unrecognised names are passed then an error is called.
%
% see https://www.e-ifu.com/search-document-metadata/OCTARAY
%
% Author: Nick Linton (2023) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications - 

% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------


if ischar(names)
    names = {names};
end

conData = getCartoConnectorElectrodeNaming();

if numel(names)==1 % maybe the user gave a Connector name
    for i = 1:numel(conData)
        if matches(conData(i).connector,names{1},IgnoreCase=true)
            namesRead = conData(i).electrodeNames;
            electrodePositions = local_readallpositions(pointFileName, conData(i));
            return %% return to user
        end
    end
end

% so now we assume that everything in names is either a bipole or a unipole
conIndex = zeros(numel(names),1);
elecIndex = zeros(numel(names),1);
bipIndex = zeros(numel(names),2);
for iName = 1:numel(names)
    for i = 1:numel(conData)
        testUNI = regexpi(names{iName},[conData(i).unipolarNaming '\d']);
        if ~isempty(testUNI) && testUNI(1)==1
            %it's a unipolar channel
            conIndex(iName) = i;
            elecIndex(iName) = local_getIndexFirstMatch(conData(i).electrodeNames, names{iName});
            if elecIndex(iName)==0
                error('Unable to match electrode name, note that multiple connector names are not allowed.')
            else
                break
            end
        end
        testBI = regexpi(names{iName},[conData(i).bipolarNaming '\d']);
        if ~isempty(testBI) && testBI(1)==1
            % it's a bipolar channel
            conIndex(iName) = i;
            tempIndex = local_getIndexFirstMatch(conData(i).bipoleNames, names{iName});
            if tempIndex==0
                error('Unable to match electrode name, note that multiple connector names are not allowed.')
            else
                bipIndex(iName,1) = local_getIndexFirstMatch(conData(i).electrodeNames, conData(i).bipoleElectrodeOneName(tempIndex));
                bipIndex(iName,2) = local_getIndexFirstMatch(conData(i).electrodeNames, conData(i).bipoleElectrodeTwoName(tempIndex));
                break
            end
        end
    end
end
validElecIndex = (elecIndex~=0);
validBipIndex = (bipIndex(:,1)~=0);
if any(~xor(validElecIndex,validBipIndex))
    error('No match or multiple matches for electrode name')
end

% get the Filenames that we will need
connectorFilenames = local_getConnectorFilenames(pointFileName);

electrodePositions = zeros(numel(names),3);
for iCon = 1:numel(conData)
    isThisConnector = (conIndex==iCon);
    if any(isThisConnector)
        connectorName = conData(conIndex).connector;
        idx = local_getIndexFirstMatch(connectorFilenames(:,1),connectorName);
        fileName = connectorFilenames{idx,2};
        [electrodeNumbering, xyz] = read_positions_on_annotation_v2(fileName);
        if ~isequal(electrodeNumbering(:), conData(conIndex).electrodeNumbers(:))
            error(['The electrode numbering is unexpected in: ' fileName])
        end
        % fill the unipolar electrode positions
        electrodePositions((isThisConnector & validElecIndex),:) = xyz(elecIndex(isThisConnector & validElecIndex) ,:);
        % take the mid-position for the bipolar electrode positions
        electrodePositions((isThisConnector & validBipIndex),:) = 0.5* ( ...
               xyz(bipIndex((isThisConnector & validBipIndex),1) ,:) + ...
               xyz(bipIndex((isThisConnector & validBipIndex),2) ,:)   );
    end
end
namesRead = names;
end

function connectorFilenames = local_getConnectorFilenames(pointFileName)
    pointExport2 = xmlgetnode(pointFileName, 'Positions', 2, 'first');
    positions = xml_readnode(pointExport2);
    [homeDir, ~,~] = fileparts(pointFileName);
    nConnectors = numel(positions.Connector);

    connectorFilenames = cell(nConnectors,2);
    count = 1;
    for iC = 1:nConnectors
        names = fieldnames(positions.Connector(iC).ATTRIBUTE);
        positionFile =  positions.Connector(iC).ATTRIBUTE.(names{1});
        if contains(positionFile,'Eleclectrode_positions_OnAnnotation.txt', 'IgnoreCase',true) %then it is an Eleclectrode_Positions_OnAnnotation.txt file
            connectorFilenames{count,1} = names{1};
            connectorFilenames{count,2} = fullfile(homeDir, positionFile);
            count = count+1;
        end
    end
    connectorFilenames(count:end,:) = [];
end


function idx = local_getIndexFirstMatch(nameList, name)
    idx = 0;
    for i = 1:numel(nameList)
        if matches(nameList{i},name)
            idx = i;
            return
        end
    end
end




          

    
    