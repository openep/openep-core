function setStimIndices(hB, iStim, isStimCaptured)
    
    if ~isequal(size(iStim), size(isStimCaptured))
        error('iStim and isStimCaptured must be the same size')
    elseif ~islogical(isStimCaptured)
        error('isStimCaptured must be a logical array')
    end
    
    hB.PrivateStimIndices = iStim;
    hB.PrivateIsStimCaptured = isStimCaptured;

end

