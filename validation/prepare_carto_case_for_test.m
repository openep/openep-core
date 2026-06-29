function [caseFolder, cleanupObj, info] = prepare_carto_case_for_test(cartoPath)
%PREPARE_CARTO_CASE_FOR_TEST Return an extracted CARTO folder for tests.
%
% Folders are returned unchanged. ZIP archives are extracted to /dev/shm when
% possible, otherwise to tempdir. The cleanup object removes extracted files.

cartoPath = char(cartoPath);
cleanupObj = onCleanup(@() []);
info = struct('sourcePath', cartoPath, 'wasArchive', false, ...
    'extractionRoot', '', 'message', '');

if isfolder(cartoPath)
    caseFolder = cartoPath;
    info.message = 'Using extracted CARTO folder.';
    return
end

if ~isfile(cartoPath)
    error('prepare_carto_case_for_test:MissingPath', ...
        'CARTO path does not exist: %s', cartoPath);
end

[~, ~, ext] = fileparts(cartoPath);
if ~strcmpi(ext, '.zip')
    error('prepare_carto_case_for_test:UnsupportedFile', ...
        'Expected a CARTO folder or ZIP file: %s', cartoPath);
end

archiveInfo = dir(cartoPath);
requiredBytes = max(archiveInfo.bytes * 4, archiveInfo.bytes + 1e9);
extractionBase = chooseExtractionBase(requiredBytes);
if isempty(extractionBase)
    error('prepare_carto_case_for_test:InsufficientSpace', ...
        ['Not enough temporary space to extract %s. ', ...
        'Need roughly %.2f GB free.'], cartoPath, requiredBytes / 1e9);
end

destination = tempname(extractionBase);
mkdir(destination);
cleanupObj = onCleanup(@() removeFolder(destination));

unzip(cartoPath, destination);
caseFolder = findExtractedCartoFolder(destination);
if isempty(caseFolder)
    error('prepare_carto_case_for_test:NoCartoFolder', ...
        'ZIP was extracted, but no CARTO study folder was found: %s', cartoPath);
end

info.wasArchive = true;
info.extractionRoot = destination;
info.message = sprintf('Extracted CARTO ZIP to %s.', destination);
end

function extractionBase = chooseExtractionBase(requiredBytes)
candidates = {};
if isfolder('/dev/shm')
    candidates{end+1} = '/dev/shm'; %#ok<AGROW>
end
candidates{end+1} = tempdir;

extractionBase = '';
for i = 1:numel(candidates)
    candidate = candidates{i};
    if usableSpaceBytes(candidate) >= requiredBytes
        extractionBase = candidate;
        return
    end
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
caseFolder = '';
meshFiles = dir(fullfile(rootFolder, '**', '*.mesh'));
for i = 1:numel(meshFiles)
    folderPath = meshFiles(i).folder;
    xmlFiles = dir(fullfile(folderPath, '*.xml'));
    names = {xmlFiles.name};
    hasStudyXml = any(~startsWith(names, '.') & ...
        ~contains(names, 'Point_Export') & ~contains(names, 'Points_Export'));
    if hasStudyXml
        caseFolder = folderPath;
        return
    end
end
end

function removeFolder(folderPath)
if isfolder(folderPath)
    rmdir(folderPath, 's');
end
end
