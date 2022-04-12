function hB = loadobj(hB)
    if isempty(hB.ChDataFileMap)
        % do nothing - we will try to reload the BardFile if egm is called
    else
        data = hB.ChDataFileMap;
        validateattributes(data, {'int16'}, {'size',[hB.NSamples hB.NChannels]});
        hB.ChDataFileMap = [];
        createEmptyMemMap(hB);
        hB.ChDataFileMap.Data.a2d = data;
    end
end

