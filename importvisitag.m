function userdata = importvisitag(userdata, dirName)
% IMPORTVISITAG provides a data structure from Carto visitag files.
%
% Usage
%   visitag = importvisitag(userdata, dirName)
%   visitag = imporvisitag()
% Where:
%   dirName is the directory with all of the files corresponding to WiseTag
%   visitag is a single data structure
%
% IMPORTVISITAG parses the data contained in a Visitag export from the
% Carto3 mapping system. The data is stored in a field '.rfindex' and the
% new userdata data structure with the appended ablation data is returned.
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
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
%
% See also, getAblationArea.m, plotAblationArea.m, plotVisitags.m

% TODO: GUI for picking the directory
% TODO: CLI for picking the directory directly from the
% userdata.<foldername> field

% read the data files
[ablSitesData, ablSitesHeader] = read_visitag_file_v1([dirName filesep() 'AblationSites.txt']);
[datPos, hdPos] = read_visitag_file_v1([dirName filesep() 'PositionsData.txt']);
[contactForceData, contactForceHeader] = read_visitag_file_v1([dirName filesep() 'ContactForceData.txt']);
[datSite, hdSite] = read_visitag_file_v1([dirName filesep() 'Sites.txt']);
[datGrid, hdGrid] = read_visitag_file_v1([dirName filesep() 'AllPositionInGrids.txt']);
[datGPos, hdGPos] = read_visitag_file_v1([dirName filesep() 'Grids.txt']);

% Access tag data
tag.X = datSite(:,lCol(hdSite, {'X' 'Y' 'Z'}));  % cartesian co-ordinates
tag.time = datSite(:,lCol(hdSite, 'DurationTime')); % seconds
tag.avgForce = datSite(:,lCol(hdSite, 'AverageForce')); % grams
tag.maxTemp = datSite(:,lCol(hdSite, 'MaxTemperature')); % C
tag.maxPower = datSite(:,lCol(hdSite, 'MaxPower')); % W
tag.Impedance.baseImp = datSite(:, lCol(hdSite, 'BaseImpedance')); %ohms
tag.Impedance.impDrop = datSite(:, lCol(hdSite, 'ImpedanceDrop')); %ohms
tag.fti = datSite(:,lCol(hdSite, 'FTI'));
tag.index.name = 'AblationIndex';
tag.index.value = datSite(:,lCol(hdSite, 'RFIndex'));

% Access grid data
% Find the unique SiteIndex (column in AllPositionInGrids)
siteIndices = datGrid(:,lCol(hdGrid, 'SiteIndex'));
uSiteIndices = unique(siteIndices);
for iSite = 1:numel(uSiteIndices)
    
    iIndex = find(siteIndices==uSiteIndices(iSite));
    tStamp = datGrid(iIndex,lCol(hdGrid, 'PosTimeStamp'));
    uniqID = datGrid(iIndex,lCol(hdGrid, 'UniqID'));
    
    uUniqID = unique(uniqID); % for example there might be 41 unique IDs for site 1
    
    for iID = 1:numel(uUniqID)
        % we now need to work with every unique uniqID
        
        % for each unique uniqID, first get the list of time stamps
        thisID = uUniqID(iID);
        data_times = tStamp(uniqID==thisID);
        
        % now use the list of timestamps to get data from the other files
        % first, locate these time stamps in PositionData
        [~, thisInd, ~] = intersect(datPos(:,lCol(hdPos, 'TimeStamp')), data_times);
        % now use the UniqID to find the X,Y,Z co-ordinates
        %identify correct row                       %identify correct column
        grid{iSite}(iID).X = datGPos( (datGPos(:,lCol(hdGPos,'UniqID'))==thisID), (lCol(hdGPos, {'X' 'Y' 'Z'})) );
        
        % now save data from PositionData
        grid{iSite}(iID).time = datPos(thisInd,lCol(hdPos,'TimeStamp'));
        grid{iSite}(iID).index.name = 'AblationIndex';
        grid{iSite}(iID).index.value = datPos(thisInd,lCol(hdPos,'RFIndex'));
        grid{iSite}(iID).impedance = datPos(thisInd,lCol(hdPos,'Impedance'));
        grid{iSite}(iID).temperature = datPos(thisInd,lCol(hdPos,'Temperature'));
        grid{iSite}(iID).power = datPos(thisInd,lCol(hdPos,'Power'));
        grid{iSite}(iID).impedanceDrop = datPos(thisInd,lCol(hdPos,'ImpedanceDrop'));
        grid{iSite}(iID).force = datPos(thisInd,lCol(hdPos,'Force'));
    end
end
userdata.rfindex.tag = tag;
userdata.rfindex.grid = grid;

% lookup column names function
    function colId = lCol(allNames, requiredName) %l for local
        if iscell(requiredName)
            for iName = 1:numel(requiredName)
                colId(iName) = find(strcmpi(allNames, requiredName{iName}));
            end
        else
            colId = find(strcmpi(allNames, requiredName));
        end
    end

end


