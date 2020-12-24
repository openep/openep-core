function [electrogramname_bip, electrogramname_uni, point_xyz_2 ] = getpointelectrogramname(point_xyz, pointFileName)
% GETELCTRODENAME finds the electrode at xyz.
%
% Usage:
%   [ electrogramname_bip, electrogramname_uni ] = getelectrogramname( point_xyz, pointFileName )
% Where:
%   point_xyz  - 
%   pointFileName  - 
%   electrogramname_bip  - the electrode name for the bipolar electrogram at point_xyz
%   electrogramname_uni - cell array of the electrode names for the
%       unipolar electrogram at point_xyz and it's associated second electrode
%       that makes up electrogramname_bip
%   point_xyz_2 - the electrode position corresponding to the second electrode
%       that makes up the bipole at point_xyz (needed for locating unipolar
%       electrograms)
%
% GETELCTRODENAME Detailed description
%
% ** ERROR noticed 30-4-2020. This program assumes sequential number of
% ** electrodes in the _onAnnotation files; however this was not the case for
% ** a Carto 3V6 case using the Pentaray connected to 20A connector. Note
% ** the numbering goes 1,2 then 1,2,3,4 repeatedly in column 1. This means
% ** that Lasso 1,2 or Lasso 3,4 were always identified as the recording
% ** electrodes, in contradiction to the header in the ECG_Export.txt file.
% **
% **       File 1-LA_P2_MAGNETIC_20_POLE_A_CONNECTOR_Electrode_Positions_OnAnnotation.txt
% **       Eleclectrode_Positions_2.0
% **       Electrode# Time X Y Z
% **       1 -1 -18.4089 -46.1829 81.5786
% **       2 -1 -22.1867 -53.66 86.9808
% **       1 -1 -13.6522 -29.4837 85.2534
% **       2 -1 -14.026 -30.9731 83.9839
% **       3 -1 -15.0154 -35.7805 80.5865
% **       4 -1 -15.3354 -37.4648 79.7378
% **       1 -1 -5.21543 -33.8042 78.6199
% **       2 -1 -6.45994 -35.0289 77.845
% **       3 -1 -11.0587 -38.1983 77.0946
% **       4 -1 -12.7287 -39.1731 77.3427
% **       1 -1 -7.75467 -49.9424 73.8663
% **       2 -1 -8.91601 -48.3531 74.1426
% **       3 -1 -12.692 -43.955 75.4862
% **       4 -1 -13.9689 -42.8323 76.1126
% **       1 -1 -12.3655 -29.5709 79.5005
% **       2 -1 -13.411 -30.7094 78.4787
% **       3 -1 -15.571 -35.4638 77.1105
% **       4 -1 -16.0402 -37.3789 77.2001
% **       1 -1 -25.5615 -30.9667 80.2972
% **       2 -1 -24.2683 -32.4533 79.986
% **       3 -1 -20.1096 -36.6017 78.8712
% **       4 -1 -18.8193 -37.9237 78.5741
%
% Author: Nick Linton (2012) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications - 
%   Steven Williams (2016)
%       Changed to 20B_1-2(139) format from 20B_1-20B_2(139) format
%       Uncomment out line 71-76 to revert back
%   Steven Williams (2017)
%       Export unipole electrograms from Carto
%   Steven Williams (2020)
%       Added description of electrode naming error
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------


pointExport2 = xmlgetnode(pointFileName, 'Positions', 2, 'first');
positions = xml_readnode(pointExport2);
 
[homeDir, ~, ~] = fileparts(pointFileName);
minDist = Inf;

closestFile = [];
closestElectrode = [];

for iC = 1:numel(positions.Connector)
    positionFile = positions.Connector(iC).ATTRIBUTE;
    names = fieldnames(positionFile);
    positionFile =  positionFile.(names{1});
    k  = strfind(lower(positionFile), 'onannotation.txt');
    if ~isempty(k) %then it is an _OnAnnotation.txt file
        positionFile = fullfile(homeDir, positionFile);
        [iElectrode, xyz] = read_positions_on_annotation_v2(positionFile);
        [iClosest, dist] = findclosestvertex(xyz, point_xyz);
        if dist(1) < minDist
            minDist = dist;
            closestElectrode = iElectrode(iClosest);
            closestFile = positionFile;
            try
                point_xyz_2 = xyz(iClosest+1,:);
            catch
                warning('xyz out of range');
            end
        end
    end
end

fileName = closestFile;
iElec = closestElectrode;

if isempty(fileName) || isempty(iElec)
    electrogramname_bip = [];
    return
end

% Now we have to translate the filename into the electrogram name that
% gives us the correct electrogram in the ECG_Export file.

identifier = {        'CS_CONNECTOR' ...
                    , 'MAGNETIC_20_POLE_A_CONNECTOR' ...
                    , 'MAGNETIC_20_POLE_B_CONNECTOR' ...
                    , 'NAVISTAR_CONNECTOR' ...
                    };
translation = {       'CS' ...
                    , '20A_' ...
                    , '20B_' ...
                    , 'M' ...
                    };
for i = 1:numel(identifier)
    if ~isempty(strfind(fileName, identifier{i}))
        electrogramname_bip = [ translation{i} num2str(iElec) '-' translation{i} num2str(iElec+1) ];
        electrogramname_uni{1} = [ translation{i} num2str(iElec) '(']; %'(' added to make sure that only one electrode is found e.g. 20B_1 could be 20B_1(61), 20B_10(70) etc
        electrogramname_uni{2} = [ translation{i} num2str(iElec+1) '('];
        if strcmpi(electrogramname_bip, '20B_7-20B_8')
            electrogramname_bip = '20B_9-20B_8';
        end
        
        % // Comment out this section to use the 20B_1-20B_2(139) rather
        % than 20B_1-2(139) format
        if i==2 || i==3 % i.e. 20A_ or 20B_
            electrogramname_bip = [ translation{i} num2str(iElec) '-' num2str(iElec+1) ];
            electrogramname_uni{1} = [ translation{i} num2str(iElec) '('];
            electrogramname_uni{2} = [ translation{i} num2str(iElec+1) '('];
            if strcmpi(electrogramname_bip, '20B_7-8')
                electrogramname_bip = '20B_9-8';
            end
        end
        % // End comment here
        return
    end
end



          

    
    