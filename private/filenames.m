function names = filenames(varargin)
% FILENAMES Gets filenames contained in a directory.
%
% Usage:
%   names = filenames(direc)
%   names = filenames(direc, searchname)
%   names = filenames(direc, searchname , ... 'option',val)
%
% Where:
%   names - a cellstr array of filenames (foldernames not included by default)
%   direc - the directory - if not passed then current directory used
%   searchname - restricts the filnames returned. Valid formats are ...
%                   name
%                   name.ext
%                   name*.ext*
%   'option',val - pairs are as follows
%       'fullfile' - true {false} - returns full file specification
%       'searchsubdirectories' - true {false} - includes all files from
%           subdirectories.
%
% FILENAMES uses DIR and then puts all the names with isdir=0 into a
% cellstr array.
%
% ALSO SEE: foldernames

% Author: Nick Linton (2011)
% Modifications - searchname changed from extension only 2013 NL
%               - 'option' added 2013 NL

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------

% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

    
    defaultDir = cd();
    defaultSearchName = '*.*';
    defaultFullFile = false;
    defaultSearchSubdirectories = false;
    
    p = inputParser;
    addOptional(p,'dir',defaultDir, @isdir);
    addOptional(p,'searchname',defaultSearchName, @ischar);
    addParameter(p,'fullfile',defaultFullFile, @islogical);
    addParameter(p,'searchsubdirectories',defaultSearchSubdirectories, @islogical);
    parse(p,varargin{:});

    
    d = dir(p.Results.dir);
    names = {d.name}';
    
    isdirec = [d.isdir];
    subdirec = names(isdirec);
    %remove '.' and '..'
    subdirec(ismember(subdirec,{'.','..'})) = [];

    names = names(~isdirec);
    
    isGood = local_checknames(names, p.Results.searchname);
    names(~isGood) = [];
    
    if p.Results.fullfile
        for i = 1:numel(names)
            names{i} = fullfile(p.Results.dir, names{i});
        end
    end
    
    if p.Results.searchsubdirectories
        for i = 1:numel(subdirec)
            subdir = fullfile(p.Results.dir, subdirec{i});
            names = vertcat( names{:} , filenames( subdir, p.Results.searchname, 'searchsubdirectories',true , 'fullfile',p.Results.fullfile  )   );
        end
    end
end



function isGood = local_checknames(names, searchName)
    [sD,sF,sE]=fileparts(searchName); %s for search
    if ~isempty(sD)
        error('The searchname includes a directory and is invalid.')
    end
    
    isGood = false(numel(names),1);
    for i = 1:numel(names)
        [~,nF,nE] = fileparts(names{i});
        if local_match(sF,nF) && local_match(sE,nE)
            isGood(i) = true;
        end
    end
end

function isMatch = local_match(s,n)
    if s(1)=='*' || strcmp(s,'.*')
        isMatch = true;
        return
    else
        iStar = find((s=='*'),1,'first');
        if isempty(iStar)
            isMatch = strcmp(s,n);
        else
            isMatch = strstartcmp(s(1:(iStar-1)), n);
        end
    end
end