function [caseFolder, cleanupObj, info] = prepare_carto_case(cartoPath)
%PREPARE_CARTO_CASE Return an extracted CARTO folder ready for import.
%
% Folders are returned unchanged. ZIP archives are extracted to /dev/shm when
% possible, otherwise to tempdir. Keep the returned cleanup object alive for
% as long as the extracted files are needed.

cartoPath = char(cartoPath);
cleanupObj = onCleanup(@() []);
info = emptyInfo(cartoPath);

if isfolder(cartoPath)
    caseFolder = cartoPath;
    info.message = 'Using extracted CARTO folder.';
    return
end

if ~isfile(cartoPath)
    error('prepare_carto_case:MissingPath', ...
        'CARTO path does not exist: %s', cartoPath);
end

[~, ~, ext] = fileparts(cartoPath);
if ~strcmpi(ext, '.zip')
    error('prepare_carto_case:UnsupportedFile', ...
        'Expected a CARTO folder or ZIP file: %s', cartoPath);
end

archiveInfo = inspect_carto_zip(cartoPath);
safetyBytes = max(0.15 * archiveInfo.uncompressedBytes, 2 * 1024^3);
requiredBytes = archiveInfo.uncompressedBytes + safetyBytes;
[extractionBase, availableBytes] = chooseExtractionBase(requiredBytes);
if isempty(extractionBase)
    error('prepare_carto_case:InsufficientSpace', ...
        ['Not enough temporary space to extract %s. Need %.2f GiB ', ...
        '(%.2f GiB archive plus safety margin).'], cartoPath, ...
        requiredBytes / 1024^3, archiveInfo.uncompressedBytes / 1024^3);
end

destination = tempname(extractionBase);
mkdir(destination);
cleanupObj = onCleanup(@() removeFolder(destination));

extractStart = tic;
extractZipArchive(cartoPath, destination);
extractionSeconds = toc(extractStart);
caseFolder = findExtractedCartoFolder(destination);
if isempty(caseFolder)
    error('prepare_carto_case:NoCartoFolder', ...
        'ZIP was extracted, but no CARTO study folder was found: %s', cartoPath);
end

info.wasArchive = true;
info.extractionRoot = destination;
info.archiveCompressedBytes = archiveInfo.compressedBytes;
info.archiveUncompressedBytes = archiveInfo.uncompressedBytes;
info.archiveFileCount = archiveInfo.fileCount;
info.requiredBytes = requiredBytes;
info.availableBytesBeforeExtraction = availableBytes;
info.extractionSeconds = extractionSeconds;
info.message = sprintf('Extracted CARTO ZIP to %s.', destination);
end

function extractZipArchive(zipPath, destination)
if isunix && isfile('/usr/bin/unzip')
    extractionLog = fullfile(destination, 'openep_unzip.log');
    command = java.util.ArrayList();
    command.add('/usr/bin/unzip');
    command.add('-q');
    command.add(zipPath);
    command.add('-d');
    command.add(destination);
    processBuilder = java.lang.ProcessBuilder(command);
    processBuilder.redirectErrorStream(true);
    processBuilder.redirectOutput(java.io.File(extractionLog));
    process = processBuilder.start();
    status = process.waitFor();
    if status ~= 0
        details = '';
        if isfile(extractionLog)
            details = strtrim(fileread(extractionLog));
        end
        error('prepare_carto_case:ExtractionFailed', ...
            'unzip failed with status %d while extracting %s. %s', ...
            status, zipPath, details);
    end
    if isfile(extractionLog)
        delete(extractionLog);
    end
else
    unzip(zipPath, destination);
end
end

function info = emptyInfo(sourcePath)
info = struct( ...
    'sourcePath', sourcePath, ...
    'wasArchive', false, ...
    'extractionRoot', '', ...
    'archiveCompressedBytes', 0, ...
    'archiveUncompressedBytes', 0, ...
    'archiveFileCount', 0, ...
    'requiredBytes', 0, ...
    'availableBytesBeforeExtraction', 0, ...
    'extractionSeconds', 0, ...
    'message', '');
end

function [extractionBase, availableBytes] = chooseExtractionBase(requiredBytes)
candidates = {};
if isfolder('/dev/shm')
    candidates{end+1} = '/dev/shm';
end
candidates{end+1} = tempdir;

extractionBase = '';
availableBytes = 0;
for i = 1:numel(candidates)
    candidate = candidates{i};
    candidateBytes = usableSpaceBytes(candidate);
    hasMemory = true;
    if strcmp(candidate, '/dev/shm')
        memoryReserveBytes = 8 * 1024^3;
        hasMemory = availableMemoryBytes() >= requiredBytes + memoryReserveBytes;
    end
    if candidateBytes >= requiredBytes && hasMemory
        extractionBase = candidate;
        availableBytes = candidateBytes;
        return
    end
end
end

function bytes = availableMemoryBytes()
bytes = Inf;
if ~isfile('/proc/meminfo')
    return
end

text = fileread('/proc/meminfo');
token = regexp(text, 'MemAvailable:\s+(\d+)\s+kB', 'tokens', 'once');
if ~isempty(token)
    bytes = str2double(token{1}) * 1024;
end
end

function bytes = usableSpaceBytes(folderPath)
bytes = 0;
try
    fileObj = java.io.File(folderPath);
    bytes = double(fileObj.getUsableSpace());
catch
    [status, out] = system(sprintf('df -Pk "%s"', folderPath));
    if status == 0
        lines = regexp(strtrim(out), '\n', 'split');
        if numel(lines) >= 2
            parts = regexp(strtrim(lines{2}), '\s+', 'split');
            if numel(parts) >= 4
                bytes = str2double(parts{4}) * 1024;
            end
        end
    end
end
end

function caseFolder = findExtractedCartoFolder(rootFolder)
candidateFolders = {};
meshFiles = dir(fullfile(rootFolder, '**', '*.mesh'));
for i = 1:numel(meshFiles)
    folderPath = meshFiles(i).folder;
    xmlFiles = dir(fullfile(folderPath, '*.xml'));
    names = {xmlFiles.name};
    hasStudyXml = any(~startsWith(names, '.') & ...
        ~contains(names, 'Point_Export') & ~contains(names, 'Points_Export'));
    if hasStudyXml
        candidateFolders{end+1} = folderPath; %#ok<AGROW>
    end
end
candidateFolders = unique(candidateFolders, 'stable');

if isempty(candidateFolders)
    caseFolder = '';
elseif isscalar(candidateFolders)
    caseFolder = candidateFolders{1};
else
    error('prepare_carto_case:AmbiguousArchive', ...
        'CARTO ZIP contains multiple study folders: %s', ...
        strjoin(candidateFolders, ', '));
end
end

function removeFolder(folderPath)
if isfolder(folderPath)
    rmdir(folderPath, 's');
end
end
