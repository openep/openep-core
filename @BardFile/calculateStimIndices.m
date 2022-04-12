function calculateStimIndices(hB)
    
    if isempty(hB.ChStim)
        hB.PrivateStimIndices = [];
        return
    end
    
    if hB.PrivateIsStimIndicesCalculated
       return
    end
    
    egm = hB.egm(':', hB.ChStim );
    
    gra = zeros(size(egm));
    gra(2:end) = (egm(2:end) - egm(1:(end-1))) * hB.SampleRate;
    
    pos = gra > hB.STIMTHRESHOLD;
    neg = gra < (-1) * hB.STIMTHRESHOLD;
    
    pos = pos(2:end) & ~pos(1:(end-1)); pos = [0 ; pos(:)];
    neg = neg(2:end) & ~neg(1:(end-1)); neg = [0 ; neg(:)];
    
    iStim = find(pos);
    
    t150 = round(hB.SampleRate * 150/1000);
    t20 = round(hB.SampleRate * 20/1000);

    %avoid any stims within 150ms of the end of the trace
    iStim(iStim>(numel(egm)-t150)) = [];
    
    for i = 1:length(iStim)
        if pos(iStim(i)) && any(neg(iStim(i):(iStim(i)+t20)))
            %then its a good bet
            pos(iStim(i) + (1:t150)) = 0;
        else
            iStim(i) = NaN;
        end
    end
    
    iStim(isnan(iStim)) = [];
    
    hB.PrivateIsStimIndicesCalculated = true;
    hB.PrivateStimIndices = iStim;
    hB.PrivateIsStimCaptured = true(size(iStim));
    verifyStimIndices(hB)

end

