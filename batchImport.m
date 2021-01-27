% ----------------------------------------------------------------------- %
% OPENEP/batchImport is a template script to automate the import of 
% datasets into OpenEP data format. This batchImport script requires a
% number of assumptions to be true in order to function.
%
%   1. Cases should be named using a study number (e.g. 001) and compressed
%   into folders named that study number.
%   2. An Excell information spreadsheet should be provided which 
%   identifies the map to import and should contain two columns where the 
%   first contains the study number and the second the number of points per 
%   map.
%
% As currently configured, the script assumes that datasets are stored on
% an external hard drive, but copies dataset zip files to the local machine
% hard drive prior to unzipping, which is usually fastest. 
%
% After importing into OpenEP dataformat, intermediate files which can
% number in the 1000s and occupy Gb of space, are removed.
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
%                             Configuration
% testmode copies and unzips files but does not process them
testmode = false; 

% Set up working directories and paths
zip_dir = '/Volumes/Extreme SSD/openep/openep_carto_export_zips';
working_dir = '/Users/Steven/Desktop/openep_working_dir';
info_spreadsheet = '/Volumes/Extreme SSD/openep/openep_datasheet.xlsx';

% Cell array of possible start of names of the study XML file
startStrings = {'CARTO' 'PVI' 'Study' '1008' 'PERS' 'AF' 'LA' 'PAF' 'AIAT'};
% ----------------------------------------------------------------------- %

% Load the study dataset Excel file
info = xlsread(info_spreadsheet, 'A:B');

% Get a list of files in zip_dir
disp('Getting list of filenames')
allFiles = nameFiles(zip_dir, 'showhiddenfiles', false);

% Loop through all the zip files
for i = 1:numel(allFiles)
    disp(['Processing case: ' allFiles{i}]);
    
    % Copy the zip file to the hard drive
    sourceFile = [zip_dir filesep() allFiles{i}];
    destinationFile = [working_dir filesep() allFiles{i}];
    status = copyfile(sourceFile, destinationFile);
    
    % Get the study number
    studyNumber = str2double(allFiles{i}(1:3));
    
    % Unzip the zip file to a folder named studyNumber in working_dir
    caseDir = ([working_dir filesep() num2str(studyNumber)]);
    unzip(destinationFile, caseDir);
    
    % Locate the carto export xml file
    allXmlFiles = nameFiles(caseDir, 'extension', 'xml', 'showhiddenfiles', false);
    iF = NaN(size(startStrings));
    for j = 1:numel(startStrings)
       temp = find(strstartcmpi(startStrings{j}, allXmlFiles)); 
       if ~isempty(temp)
           iF(j) = temp;
       end
    end
    
    if nnz(~isnan(iF)) > 1
        % Multiple possible files found so the user has to choose which xml file to load
        iF(isnan(iF)) = [];
        beep()
        disp([num2str(numel(iF)) ' possible XML files found:'])
        for j = 1:numel(iF)
           disp([ num2str(j) '    ' allXmlFiles{iF(j)}]);
        end  
        result = input('Enter choice: ');
        studyXmlFile = [case_dir filesep() allXmlFiles{iF(result)}];
    else % Automatically select the relevant xml file
        iF(isnan(iF)) = [];
        studyXmlFile = [caseDir filesep() allXmlFiles{iF}];
    end
    
    % Check that an XML file has been identified
    if ~isfile(studyXmlFile)
       error(['RUN_PEPR_EXPERIMENT: Specified XML file - ' studyXmlFile ' - does not exist']); 
    end
    
    % Read the number of points from the spreadsheet
    numpts = info(info(:,1)==studyNumber,2);
    
    % Specify the output directory and filename
    outputFile = ([working_dir filesep() [num2str(studyNumber) '.mat']]);
    
    % Start a timer
    tic();
    
    % Run importcarto_mem from the command line
    if ~testmode
        userdata = importcarto_mem(studyXmlFile ...
            , 'maptoread', numpts ...
            , 'refchannel', 'CS9-CS10' ...
            , 'ecgchannel', 'V1' ...
            , 'savefilename', outputFile ...
            );
    end
    
    % Stop a timer
    elapsedTime = toc();
    
    % Get the size of caseDir and destinationFile
    sizeUnZipped = du(['-sh ' caseDir]);
    sizeZipped   = du(['-sh ' destinationFile]);
    
    % Remove the carto files and zip file
    deleteFolder(caseDir);
    delete(destinationFile);
    
    % Dump the time and filesizes to a text file
    fid = fopen([working_dir filesep() [num2str(studyNumber) '.data']], 'wt');
    fprintf(fid, [num2str(elapsedTime) '\n']);
    fprintf(fid, [sizeUnZipped '\n']);
    fprintf(fid, [sizeZipped '\n']);
    fclose(fid);
    
end