function createEmptyMemMap(hB)
% creates an empty memory map and assigns it

    if ~isempty(hB.ChDataFileMap)
        error('@BARDFILE/CREATEEMPTYMEMMAP: There already exists a memmapfile.')
    end


    %calculate size of data
    dataSize = hB.NChannels * hB.NSamples;

    %create a binary file large enough to hold all of the data
    dataFileName = tempname;

    [dataFileID, message] = fopen( dataFileName, 'w+' );
    if dataFileID == (-1)
        disp('BARDFILE.m: There was a problem creating a temporary file:')
        error(message)
    end
    disp(['A temporary file was created: ' dataFileName])

    fwrite(dataFileID, 0, 'int16', (dataSize-1)*2); %each number is 2 bytes
    fclose(dataFileID);

    %now create a memory map to the file
    hB.ChDataFileMap = memmapfile( dataFileName, ...
        'Format', { 'int16' [hB.NSamples hB.NChannels] 'a2d' }, ...
        'Writable', true ...
        ); %'a2d to signify it is the raw output from the A/D converter.

end