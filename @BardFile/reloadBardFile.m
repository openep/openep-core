function reloadBardFile(hB)
% @BARDFILE/reloadBardFile reloads the egm data from file into hB.ChDataFileMap.
% Usage:
%   reloadBardFile(hB)

% Author: Nick Linton (2013)
% Modifications - 

    persistent masterDir
    masterDir = getApsafBardRawDataPath;
    
    if ~isempty(hB.ChDataFileMap)
        error('@BARDFILE/LOADBARDFILE: there is already a datafilemap - write code to delete it if you want to do this.')
    end
    
    % make a temporary copy of the BardFile (or its suubclass object)
    hBNew = clone(hB);
    
    % check that the filename is ok - if not then search for it
    if ~isfile(hBNew.FileName)
        n = hBNew.FileName;
        n(n=='\') = filesep();
        [~,n,e] = fileparts(n);
        disp(masterDir);
        if ~isdir(masterDir)
            h = msgbox([n e ' could not be found. Press OK to choose a master directory. The file will be searched for in this and subdirectories.'], 'CreateMode', 'modal');
            uiwait(h);
            masterDir = uigetdir(cd(),'Choose master directory for future BardFiles.');
        end
        fnames = filenames(masterDir, [n e], 'fullfile',true , 'searchsubdirectories',true );
        if numel(fnames)==1
            warning([n e ' taken from the Master Directory.'])
            hBNew.FileName = fnames{1};
        else
            error([n e ' not found'])
        end
    end
    
    % load the Bardfile
    loadBardFile(hBNew);
    
    % make sure nothing has changed ...
    hBNew.FileName = hB.FileName;
    meta = ?BardFile;
    ndp = findobj([meta.PropertyList],'Dependent',false);
    
    for idx = 1:length(ndp)
        name = ndp(idx).Name;
        if ~isequal(hBNew.(name), hB.(name))
            % we have a 'problem'
            if ~strcmpi(name, 'ChDataFileMap')
                disp(':')
                disp(':')
                disp(':')
                disp(':')
                disp('BARDFILE/reloadBardFile: When reloading the file')
                disp(['the value of BardFileObject.' name ' changed. The file does not'] )
                disp('seem compatible with the BardFile object. No values have been changed.' )
                disp('Use "loadBardFile" instead if this is what you intended.' )
                error(['Problem loading: ' hB.ShortFileName '  See message printed above.'])
            end
        end
    end


    
    % so all properties are ok
    hB.ChDataFileMap = hBNew.ChDataFileMap;
    
    % now get rid of hBNew without deleting the memmap file
    hBNew.ChDataFileMap = [];
    clear hBNew

end


