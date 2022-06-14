function nU = openep_updateuserdata(f, varargin)
% OPENEP_UPDATEUSERDATA is used to modify userdata
%
% Usage:
%   newUserData = openep_updateuserdata(dir, varargin)
%
% Where:
%   f      - path to a single userdata .mat file, or userdata
%   nU     - the new userdata
%
% OPENEP_UPDATEUSERDATA accepts the following parameter-value pairs
%   'outputpath'    {[]} | <string>
%   'version'       {'-v7'} | '-v7.3' | '-v6'
%   'unpacktrirep'  {false} | true
%
% OPENEP_UPDATEUSERDATA is used to make non-data related changes to
% userdata. The following interventions are possible
% 1) By specifying 'version', the save version can be modified (see help 
%    save). Note that unless otherwise specified, '-v7.3' is used, and that
%    saving to disc is optional.
% 2) Convert userdata.surface.triRep to a structure rather than a triRep.
%    This is useful for example if userdata is subsequently to be loaded
%    into a non-Matlab environment
% Note that this function handles single files only. Paths to single files 
% are passed in as the argument, f. Userdata already loaded can also be
% passed in as the argument, f. Limited error check is done on f.
%
% Author: Steven Williams (2021)
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

nStandardArgs = 1; % UPDATE VALUE
versionString = '-v7';
unPackTriRep = false;
outputPath = [];
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'version'
                versionString = varargin{i+1};
            case 'unpacktrirep'
                unPackTriRep = varargin{i+1};
            case 'outputpath'
                outputPath = varargin{i+1};
        end
    end
end

% Load or access the userdata
try
    if ischar(f)
        if isfile(f)
            load(f,'userdata');
        end
    elseif isstruct(f)
        userdata = f;
    else
        error('OPENEP/OPENEP_UPDATEUSERDATA: Invalid input, no userdata found')
    end
catch
    error('OPENEP/OPENEP_UPDATEUSERDATA: Invalid input, no userdata found')
end

% Unpack the TriRep (for example, for loading the .mat file into Python
if unPackTriRep
    tr = getMesh(userdata);
    X = tr.X;
    Triangulation = tr.Triangulation;
    userdata.surface.triRep = [];
    userdata.surface.triRep.X = X;
    userdata.surface.triRep.Triangulation = Triangulation;
    
    % Store comment about triRep
    userdata.notes{end+1} = [date ': TriRep unpacked using openep_updateuserdata.m'];
end

if ~isempty(outputPath)
    % Store comment about file format
    userdata.notes{end+1} = [date ': Converted to ' versionString ' using openep_updateuserdata.m'];
    save(outputPath, 'userdata', versionString);
    disp(['OPENEP/OPENEP_UPDATEUSERDATA: Data saved to ' outputPath]);
end

% output the updated userdata to the calling environment
nU = userdata;

end
