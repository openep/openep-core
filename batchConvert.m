function batchConvert(inputDir, outputDir, varargin)
% batchConvert Converts original style OpenEP datasets into the new format.
%
% In the newer format, userdata.surface.triRep is a regular struct rather
% than a TriRep. This is to allow for interopability with OpenEP-Py, which
% cannot load TriRep objects. This function also allows options to
% pre-compute normals.
%
% Usage:
%   batchConvert(inputDir, outputDir)
%
% Where:
%   inputDir - the full path to the directory containing the .mat files to be
%              converted.
%   outputDir - the directory to which the new files will be written. The filenames
%               will be the same as those in inputDir.
%
% Author: Paul Smith (2021) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%       Steven Williams (2022) Updated documentation and added notes
%       Steven Williams (2023) Added option to pre-compute normals
%
% Info on Code Testing:
% ---------------------------------------------------------------
% batchConvert(inputDir, outputDir, 'computenormals', true)
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% parse input arguments
nStandardArgs = 2;
computeNormals = false;
if nargin > nStandardArgs
    for i = nStandardArgs:2:nargin-1-nStandardArgs
        switch lower(varargin{i})
            case 'computenormals'
                computeNormals = varargin{i+1};
        end
    end
end

% Create the output folder if it does not exist
if ~exist(outputDir, 'dir')
   mkdir(outputDir)
end
 
% Get a list of files
disp('Getting list of filenames')
allFiles = nameFiles(inputDir, 'showhiddenfiles', false, 'extension', 'mat');

% Iterate through each file and convert the OpenEP dataset to the newer
% format
for i = 1:numel(allFiles)

   disp(['working on file:  ' allFiles{i}])
   load([inputDir filesep() allFiles{i}]);

   t.X = userdata.surface.triRep.X;
   t.Triangulation = userdata.surface.triRep.Triangulation;
   userdata.surface.triRep = t;

   % Store comment about what we have done
   userdata.notes{end+1} = [date ': data set converted using batchConvert.m'];

   if computeNormals
       disp('precomputing normals')
       [~, userdata] = getNormals(userdata);

       % Store comment about what we have done
       userdata.notes{end+1} = [date ': normals added using batchConvert.m'];
   end
   
   % We save as -v7 because it's faster to load in OpenEP-py than -v7.3,
   % and the saved file is significantly smaller compared to -v6 files.
   outputFile = [outputDir filesep() allFiles{i}];
   disp(['saving file: ' outputFile])
   save(outputFile, 'userdata', '-v7');

end

end
