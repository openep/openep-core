function hBSave = saveobj(hB)
% we won't save the egm data but leave that in its file
    hBSave = clone(hB);
    if hB.FullSave
        hBSave.ChDataFileMap = hB.ChDataFileMap.Data.a2d;
    else
        hBSave.ChDataFileMap = [];
    end
end
