function b = egm(hB, iTime, iChannel)
% @BARDFILE/egm     Accesses a BARDFILE electrogram.
% Usage:
%   b = egm(bardFileObj)
% Author: Nick Linton (2009)
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------

if isempty(hB.ChDataFileMap)
    %then try to reload it
    reloadBardFile(hB)
    if isempty(hB.ChDataFileMap)
        % if its still empty then nothing to do ...
        b = [];
        return
    end
end
                        
if ischar(iChannel)  ||  iscellstr(iChannel)
    iChannel = chNames2Indices(hB, iChannel);
end

if ischar(iTime) && strcmpi(iTime, ':')
    b = double( hB.ChDataFileMap.Data.a2d(:,iChannel) );
else
    b = double( hB.ChDataFileMap.Data.a2d(iTime,iChannel) );
end

b = b * 2 / (2^hB.ADCBITS);
for i = 1:length(iChannel)
    b(:,i) = b(:,i) * hB.ChRange(iChannel(i)) ;
end
