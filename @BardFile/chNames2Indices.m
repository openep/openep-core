function ind = chNames2Indices(hB, chNames)
% if the user supplies a list of Labels, instead of indices, then these are
% converted. For example, if the user types obj.egm({'HISd' 'HISp'}) and the labels
% in obj.cLabels are { 'I', 'aVF', 'HISd', 'HISp' } then ind = [3,4].

    if ~iscellstr(chNames) && ~ischar(chNames)
        error('BARDFILE/FINDCHANNELINDEX: unrecognised format.')
    end
    if ischar(chNames)
        if strcmpi(chNames, ':')
            ind = 1:hB.NChannels;
            return
        else
            chNames = {chNames};
        end
    end
    if iscellstr(chNames)
        if numel(chNames)==1 && strcmpi(chNames{1}, '')
            ind = [];
        else
            ind = zeros(1,length(chNames));
            for i = 1:length(chNames)
                match = strcmp(chNames{i}, hB.ChName);
                if sum(match) == 0
                    error('@BARDFILE/SUBSREF: No match was found for one of the specified labels.')
                elseif sum(match >= 2)
                    error('@BARDFILE/SUBSREF: More than one match was found for one of the specified labels.')
                end
                ind(i) = find(match);
            end
        end
    end
end