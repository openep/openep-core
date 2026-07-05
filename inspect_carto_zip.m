function info = inspect_carto_zip(zipPath)
%INSPECT_CARTO_ZIP Read archive sizes from the ZIP central directory.

zipPath = char(zipPath);
if ~isfile(zipPath)
    error('inspect_carto_zip:MissingFile', ...
        'ZIP file does not exist: %s', zipPath);
end

[~, ~, ext] = fileparts(zipPath);
if ~strcmpi(ext, '.zip')
    error('inspect_carto_zip:UnsupportedFile', ...
        'Expected a ZIP file: %s', zipPath);
end

fileInfo = dir(zipPath);
info = struct( ...
    'compressedBytes', double(fileInfo.bytes), ...
    'uncompressedBytes', 0, ...
    'fileCount', 0);

zipFile = java.util.zip.ZipFile(java.io.File(zipPath));
cleanupObj = onCleanup(@() zipFile.close());
entries = zipFile.entries();
unknownSize = false;

while entries.hasMoreElements()
    entry = entries.nextElement();
    if entry.isDirectory()
        continue
    end

    entrySize = double(entry.getSize());
    if entrySize < 0
        unknownSize = true;
    else
        info.uncompressedBytes = info.uncompressedBytes + entrySize;
    end
    info.fileCount = info.fileCount + 1;
end

if unknownSize
    error('inspect_carto_zip:UnknownEntrySize', ...
        'ZIP contains entries whose uncompressed size is unavailable: %s', zipPath);
end
end
