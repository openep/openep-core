function [electrodePositions, namesRead, electrodePositionsAll ] = read_electrodePositionsOnAnnotation(names, pointFileName)
% READ_ELECTRODEPOSITIONSONANNOTATION gets positions of named electrodes.
% Usage:
%   [electrodePositions, namesRead, electrodePositionsAll ] = ...
%                 read_electrodePositionsOnAnnotation(names, pointFileName)
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
%   electrodePositions - a double array of positions, size [numel(namesRead),3]
%   namesRead - usually the same as names but [] is returned in the case that
%               no match is found for a particular name.
%   electrodePositionsAll - double array of the positions of all electrodes
%                           on connectors shared by any electrode referred
%                           to in 'names', size [nElectrodes,3].
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
        if matches(conData(i).connector,names{1})%,IgnoreCase=true)
            namesRead = conData(i).electrodeNames;
            electrodePositions = local_readallpositions(pointFileName, conData(i));
            return %% return to user
        end
    end
end

% Assume that everything in names is either a bipole or a unipole.
% Find which connector each name is associated with and if it is bi- or
% uni-
conIndex = zeros(numel(names),1);
isBipolar = false(numel(names),1);
for iName = 1:numel(names)
    for i = 1:numel(conData)
        testUNI = regexpi(names{iName},[conData(i).unipolarNaming '\d']);
        if ~isempty(testUNI) && testUNI(1)==1
            %it's a unipolar channel
            %isBipolar(iName) = false;
            conIndex(iName) = i;
            break
        end
        testBI = regexpi(names{iName},[conData(i).bipolarNaming '\d']);
        if ~isempty(testBI) && testBI(1)==1
            % it's a bipolar channel
            isBipolar(iName) = true;
            conIndex(iName) = i;
            break
        end
    end
end



% get the Filenames that we will need
connectorFilenames = local_getConnectorFilenames(pointFileName);

electrodePositions = zeros(numel(names),3);
electrodePositionsAll = [];
for iCon = 1:numel(conData)
    isThisConnector = (conIndex==iCon);
    if any(isThisConnector)
        namesThisConnector = names(isThisConnector);
        connectorName = conData(conIndex).connector;
        
        idx = local_getIndexFirstMatch(connectorFilenames(:,1),connectorName);
        fileName = connectorFilenames{idx,2};
        [electrodeNumbering, xyz] = read_positions_on_annotation_v2(fileName);
        
        electrodePositionsAll = [electrodePositionsAll ; xyz]; %#ok<AGROW>
        
        expectedElectrodeNumbering = conData(conIndex).electrodeNumbers(:);
        expectedElectrodeNames = conData(conIndex).electrodeNames;
        
        if ~isequal(electrodeNumbering(:), expectedElectrodeNumbering)
            % try removing the optionalElectrodes
            expectedElectrodeNumbering(conData(conIndex).optionalElectrodes) = [];
            expectedElectrodeNames(conData(conIndex).optionalElectrodes) = [];
            if ~isequal(electrodeNumbering(:), expectedElectrodeNumbering)
                warning(['NAN returned: The electrode numbering cannot be resolved in: ' fileName])
                electrodePositions(isThisConnector,:) = nan(numel(namesThisConnector),3);
                names(isThisConnector) = cell(size(isThisConnector));
                continue %to the next iCon value in the loop
            end
        end
        
        isBipolarThisConnector = isBipolar(isThisConnector);
        electrodePositionsThisConnector = zeros(numel(namesThisConnector),3);
        for iName = 1:numel(namesThisConnector)
            if isBipolarThisConnector(iName)
                idx = local_getIndexFirstMatch(conData(i).bipoleNames, namesThisConnector{iName});
                if idx==0; error('Unable to match electrode name'); end
                idx1 = local_getIndexFirstMatch(expectedElectrodeNames, conData(i).bipoleElectrodeOneName(idx));
                idx2 = local_getIndexFirstMatch(expectedElectrodeNames, conData(i).bipoleElectrodeTwoName(idx));
                % take the mid-position for the bipolar electrode positions
                electrodePositionsThisConnector(iName,:) = 0.5* ( ...
                    xyz(idx1 ,:) + ...
                    xyz(idx2 ,:)   );
            else %it's unipolar
                idx = local_getIndexFirstMatch(expectedElectrodeNames, namesThisConnector{iName});
                if idx==0; error('Unable to match electrode name'); end
                electrodePositionsThisConnector(iName,:) = xyz(idx ,:);
            end
        end
        electrodePositions(isThisConnector,:) = electrodePositionsThisConnector;
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
        if contains(positionFile,'Eleclectrode_positions_OnAnnotation', 'IgnoreCase',true) %then it is an Eleclectrode_Positions_OnAnnotation.txt file
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




          

    
    