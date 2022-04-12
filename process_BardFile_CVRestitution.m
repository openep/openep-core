function DATA = process_BardFile_CVRestitution(bardFilePath) %#ok<INUSD>
% PROCESS_BARDFILE_CVRESTITUTION Creates CV restitution and ERP restitution
% data from S1S2 pacing protocols
%
% Usage:
%   DATA = process_BardFile_CVRestitution(fold)
% Where:
%   bardFilePath - is the directory to be processed, or 'openfile'
%   DATA - is the output data [S1S2 LAT]
%
% PROCESS_BARDFILE_CVRESTITUTION Detailed description goes here
%
% Author: Steven Williams (2016)
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% 
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% parse input data, get the directory
if strcmpi(bardFilePath, 'openfile')
    bardFilePath = uigetdir('/Users/', 'Select directory of BARD files ...');
    if bardFilePath==0
        disp('Operation cancelled');
        return
    end
end
files = nameFiles(bardFilePath);

% preassign the output
DATA = [NaN NaN];

% create a figure
figure

% process the BARD files
for iFile = 2:numel(files)
    
    % load the BARD file
    hB = BardFile([bardFilePath filesep files{iFile}]);
    
    % first time only, choose the output channel
    if iFile==2
        [selection,ok] = listdlg('ListString', hB.ChName ...
            , 'SelectionMode', 'single' ...
            , 'ListSize', [300, 600] ...
            , 'InitialValue', 1 ...
            , 'PromptString', 'Select the recording channel'  ...
            );
        if ~ok
            disp('Operation cancelled')
            return
        end
        electrode = hB.ChName{selection};
    end
    
    % Work out the S1S2 coupling interval from the filename
    S1S2 = str2double(files{iFile}(5:8));
    
    % Get the extrastimulus beat
    S1 = hB.beat(hB.NStim,electrode);
    
    % Cancel the stim artefact
    s1 = S1; s1(1:16) = 0;
    
    % Apply the non-linear energy operator
    [nl, nlfilt, act] = nleo(s1); %#ok<ASGLU>
    
    % Plot the calculated activation times for checking
    if iFile>2
        set([egmPlot actPlot], 'color', [.8 .8 .8]);
    end
    egmPlot = plot(S1, 'k'); hold on;
    actPlot = plot(act, 'r');
    title([electrode ' S1S2=' num2str(S1S2)]);
    xlabel('LAT (ms)');
    ylabel('Bipolar Voltage (mV)');
    if iFile==2
        legend({'EGM', 'Activation'})
    end
    set(gca, 'ylim', [-2 2]);
    
    % save the data
    lat = find(act);
    disp(['S1S2=' num2str(S1S2) 'ms, activation time ' num2str(lat(1)) 'ms']);
    DATA(iFile,:) = [S1S2 lat(1)];
    
    % wait user input
    pause
end

end