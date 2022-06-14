function batchAnonymise(inputDir, outputDir)
% batchAnonymise Anonymises userdata.
%
% In order to anonymise userdata, the fields containing the path to the
% original data are removed as we cannot guarantee this do not contain
% potentially identifiable data in other centres.
%
% Usage:
%   batchAnonymise(inputDir, outputDir)
%
% Where:
%   inputDir - the full path to the directory containing the .mat files to be
%              anonymised.
%   outputDir - the directory to which the new files will be written. The filenames
%               will be the same as those in inputDir.
%
% Author: Steven Williams (2021) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%

% Create the output folder if it does not exist
if ~exist(outputDir, 'dir')
   mkdir(outputDir)
end

% Get a list of files
disp('Getting list of filenames')
allFiles = nameFiles(inputDir, 'showhiddenfiles', false, 'extension', 'mat');

% Iterate through each file and anonymise the OpenEP dataset
for i = 1:numel(allFiles)

   disp(['working on file:  ' allFiles{i}])
   load([inputDir filesep() allFiles{i}]);

   if isfield(userdata, 'cartoFolder')
       userdata.cartoFolder = [];
   end

   % Store comment about what we have done
   if ~isfield(userdata, 'notes')
       userdata.notes = [];
   end
   userdata.notes{end+1} = [date ': data set anonymised using batchAnonymise.m'];

   % We save as -v7 because it's faster to load in OpenEP-py than -v7.3,
   % and the saved file is significantly smaller compared to -v6 files.
   outputFile = [outputDir filesep() allFiles{i}];
   disp(['saving file: ' outputFile])
   save(outputFile, 'userdata', '-v7');

end

end
