function deleteFolder( folderPath, varargin )
% DELETEFOLDER Deletes the folder folderPath
%
% Usage:
%   deleteFolder( folderPath )
% Where:
%   folderPath - the full path to the folder
%
% DELETEFOLDER Detailed description goes here
%
% DELETEFOLDER accepts the following parameter-value pairs
%   'verbose'     {false}|true
%
% Author: Steven Williams (2020)
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

% Parse input
nStandardArgs = 1; % UPDATE VALUE
verbose = false;
if nargin > nStandardArgs
    for i = nStandardArgs+1:2:nargin
        switch varargin{i}
            case 'verbose'
                verbose = varargin{i+1};
        end
    end
end

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(folderPath)
  error('DELETEFOLDER: Specified path does not exist')
  return;
end

% Get a list of all files in the folder with the desired file name pattern.
theFiles = nameFiles(folderPath);

% Delete each file one by one
for k = 1 : length(theFiles)
  fullFileName = fullfile(folderPath, theFiles{k});
  if verbose
    fprintf(1, 'Now deleting %s\n', fullFileName);
  end
  delete(fullFileName);
end

% Remove the directory
rmdir(folderPath)

end