% ----------------------------------------------------------------------- %
% OPENEP/batchConvert is a template script which can be used to convert
% original style OpenEP datasets into the new format.
%
% In the newer format, userdata.surface.triRep is a regular struct rather
% than a TriRep. This is to allow for interopability with OpenEP-Py, which
% cannot load TriRep objects.
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
%                             Configuration
% Select the full path to the directory containing the .mat files to be
% converted.
working_dir = '/Users/paul/Documents/openep/openep-core';

% The files will be written to a different directory. The same filenames
% will be used.
output_dir = '/Users/paul/Documents/openep/new-openep-files';
% ----------------------------------------------------------------------- %

% Create the output folder if it does not exist
if ~exist(output_dir, 'dir')
   mkdir(output_dir)
end

% Get a list of files
disp('Getting list of filenames')
allFiles = nameFiles(working_dir, 'showhiddenfiles', false, 'extension', 'mat');

% Iterate through each file and convert the OpenEP dataset to the newer
% format
for i = 1:numel(allFiles)

   disp(['working on file:  ' allFiles{i}])
   load([working_dir filesep() allFiles{i}]);

   t.X = userdata.surface.triRep.X;
   t.Triangulation = userdata.surface.triRep.Triangulation;
   userdata.surface.triRep = t;

   outputFile = [output_dir filesep() allFiles{i}];
   disp(['saving file: ' outputFile])
   save(outputFile, 'userdata', '-v6');  % v6 is faster to load in OpenEP-Py

end
