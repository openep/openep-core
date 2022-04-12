classdef BardFile  < hgsetget
% BARDFILE is a class set up to process data from the Bard.
% The Bard exports text files and this class has been created to read them.
% The electrograms are then stored on disk and a memory map to them is
% created. This allows very large files to be handled because the data does
% not need to be loaded into dynamic memory.
% The following examples indicate how a BardFile can be loaded and used:
%   x = BardFile('filename') - loads the data into a BardFile object, x
%   disp(x) (or typing 'x') - displays information about it
%   x.egm(Label) - returns the electrogram(s) corresponding to Label. If
%       Label is 'all' then all electrograms are returned in an array.
% Author: Nick Linton (2009)
% Modifications - 



% Public properties
    properties (SetAccess = 'protected') %orinally CONSTANTS but could be a nightmare in future.
        ADCBITS = 16;   %The number of bits in the A/D Converter
        STIMTHRESHOLD = 1000; %5V/ms
    end
    properties %(SetAccess = 'protected')
        FileName = [];
        StartTime = [];
        NChannels = [];        
        NSamples = [];
        SampleRate = [];
        ChName = [];    %prefix 'c' used for channel specific data
        ChRange = [];
        ChLow = [];
        ChHigh = [];
    end
    properties (Access = 'private')
        PrivateStimIndices = [];
        PrivateIsStimCaptured = [];
        PrivateIsStimIndicesCalculated = false;
        PrivateChStim = [];
        PrivateFilter = [100 500];
        PrivateDFilt = [];
    end
    properties %allow full access
        FullSave = false; %set this to true if you want to save all of the egms, false if you want the computer to reload egms from Bard text file
    end
    properties (Dependent = true)
        DFilt
        Filter
        RedoStimIndices
        ChStim
        ChStimName
        StimIndices
        StimCapturedIndices
        IsStimCaptured
        NStim
        PacedEgm
        NyquistFreq
        ShortFileName
    end
    properties (SetAccess = 'protected') %make it's GetAccess 'protected' too when debugged.
        ChDataFileMap = [];
    end


% Class methods
    methods (Static = true)
        hB = loadobj(hB)
    end
    methods %generic
        hBSave = saveobj(hB)
        newObj = clone(oldObj)
    end
    methods
        ind = chNames2Indices(hB, chNames)
        calculateStimIndices(hB)
        loadBardFile(hB)
        reloadBardFile(hB)
        [fe, e] = filtEgm(hP, iTime, iChannel, varargin)
        [e, precedingCL, startIndex] = beat(hB, nBeat, iChannel, varargin)
    end
    % Get Methods
    methods
        function nf = get.NyquistFreq(hB)
            nf = hB.SampleRate/2;
        end
        function ch = get.ChStim(hB)
            ch = hB.PrivateChStim;
        end
        function name = get.ChStimName(hB)
            ch = hB.PrivateChStim;
            if isempty(ch)
                name = '';
            else
                name = hB.ChName{ch};
            end
        end
        function name = get.ShortFileName(hB)
            fileNameString = hB.FileName;
            fileNameString(fileNameString=='\')='/';
            [~,name,ext] = fileparts(fileNameString);
            name = [name ext];
        end
        function iStim = get.StimIndices(hB)
            if isempty(hB.PrivateStimIndices)
                calculateStimIndices(hB);
            end
            iStim = hB.PrivateStimIndices;
        end
        function isCaptured = get.IsStimCaptured(hB)
            if isempty(hB.PrivateIsStimCaptured)
                calculateStimIndices(hB);
            end
            isCaptured = hB.PrivateIsStimCaptured;
        end
        function iStimCap = get.StimCapturedIndices(hB)
            iStimCap = hB.StimIndices(hB.IsStimCaptured);
        end
        function n = get.NStim(hB)
            n = numel(hB.StimIndices);
        end
        function filt = get.DFilt(hP)
            if isempty(hP.PrivateDFilt)
                hP.Filter = hP.PrivateFilter;
            end
            filt = hP.PrivateDFilt;
        end
        function f = get.Filter(hP)
            f = hP.PrivateFilter;
        end
       
    end % "get" methods
    
    % Set Methods
    methods
        function set.ChStim(hB, ch)
            if strcmpi(ch, 'none')
                hB.PrivateChStim = [];
            else
                ind = chNames2Indices(hB, ch);
                if ~isequal(ind, hB.PrivateChStim)
                    hB.PrivateChStim = ind;
                end
            end
        end
        function set.Filter(hB, bandHz)
            validateattributes(bandHz, {'numeric'}, {'numel',2});
            if ~isequal(hB.PrivateFilter, bandHz(:)') || isempty(hB.PrivateDFilt)
                hB.PrivateFilter = bandHz(:)';
                w = bandHz / hB.NyquistFreq;
                [z,p,k] = butter(2, w);
                [sos,g] = zp2sos(z,p,k);
                hB.PrivateDFilt = dfilt.df2sos(sos,g);
            end
        end
        function set.StimIndices(hB, iStim)
            validateattributes(iStim, {'numeric'}, {'vector','integer', '>=',1 , '<=',hB.NSamples});
            hB.PrivateStimIndices = iStim;
        end
        function set.IsStimCaptured(hB, isCaptured)
            validateattributes(isCaptured, {'numeric'}, {'vector','integer', '>=',0 , '<=',1 , 'size',size(hB.StimIndices)});
            hB.PrivateIsStimCaptured = logical(isCaptured);
        end
    end %"set" methods
    
    
    % BardFile method
    methods
        function delete(hB)
            if ~isempty(hB.ChDataFileMap) && isa(hB.ChDataFileMap, 'memmapfile') %the ChDataFileMap may have been set to int16 as part of saveobj
                fname = hB.ChDataFileMap.Filename;
                hB.ChDataFileMap = [];
                disp(['A temporary file was destroyed: ' fname])
                delete(fname)
            end
        end
        % locally written methods
        function obj = BardFile(varargin)
            persistent openDir
            if nargin == 0
                return
            elseif nargin == 1
                if ischar(varargin{1}) && strcmpi(varargin{1}, 'openfile')
                    if isempty(openDir)
                        openDir = matlabroot();
                    end
                    if ~isdir(openDir); openDir = matlabroot(); end
                    
                    [filename,pathname] = uigetfile('*.txt' ...
                        , 'Select the Bard export file.', openDir ...
                        ,  'MultiSelect', 'on' ...
                        );
                    if isnumeric(filename)  &&  filename == 0
                        return
                    else
                        openDir = pathname;
                    end
                    if ischar(filename)
                        allnames = {fullfile(pathname, filename)};
                    elseif iscellstr( filename )
                        n = numel(filename);
                        allnames = cell(n,1);
                        for i = 1:n
                            allnames{i} = fullfile(pathname, filename{i});
                        end
                    end
                else
                    allnames = varargin{1};
                    if ischar(allnames)
                        allnames = {allnames};
                    end
                end
                
                n = numel(allnames);
                if n>1
                    if ~isa(obj,'BardFile')
                        error('BARDFILE/BARDFILE: if a subclass calls the Bardfile constructor, then the constructor cannot change the dimension of the object.')
                    end
                    obj(n,1) = BardFile();
                end
                for iObj = 1:n
                    
                    name = allnames{iObj};
                    
                    if isempty(name)
                        return
                    end
                    [pathname, filename, ext] = fileparts(name);    %but a path may not have been specified
                    if isempty(pathname)
                        pathname = cd();
                        name = fullfile(pathname, [filename ext]);
                    end
                    
                    obj(iObj).FileName = name; %#ok<*AGROW>
                    loadBardFile(obj(iObj))
                end
            end
        end
    end % methods
%         
end % classdef

