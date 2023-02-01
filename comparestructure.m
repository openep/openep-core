function [result, commentary] = comparestructure(a,b)
% COMPARESTRUCTURE.
% Usage:
%   [result, commentary] = comparestructure(a,b)
% Inputs:
%   a,b  - input structures to compare
% Outputs:
%   result  - logical
%   commentary - a cellarray of text messages explaining the differences
%
% Author: Nick Linton (2023)
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------

% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------


    local_addmsg(0, 'reset');

    if isequal(a,b)
        result = true;
        commentary = '';
    	return
    else
        result = false;
        local_compare(a,b,0)
        commentary = local_addmsg(0, 'retrieve');
    end
end
            
function local_compare(a,b,level)
    fieldnamesA = sort(fieldnames(a));
    fieldnamesB = sort(fieldnames(b));
    fieldNamesAll = unique([fieldnamesA;fieldnamesB]);
    
    badNamesA = false(size(fieldNamesAll));
    badNamesB = false(size(fieldNamesAll));
    if ~isequal(fieldnamesA, fieldnamesB)
        for i = 1:numel(fieldNamesAll)
            if ~any(matches(fieldnamesB, fieldNamesAll{i}))
                badNamesA(i) = true;
                local_addmsg(['.' fieldnamesA{i} ' is not in the second structure']);
            end
            if ~any(matches(fieldnamesA, fieldNamesAll{i}))
                badNamesB(i) = true;
                local_addmsg(['.' fieldnamesA{i} ' is not in the second structure']);
            end
        end
    end
    
    level = level + 1;
% now go through the fieldNames that are shared.
    sharedFieldNames = fieldNamesAll(~badNamesA & ~badNamesB);
    for i = 1:numel(sharedFieldNames)
        f = sharedFieldNames{i};
        if ~strcmp(class(a.(f)),class(b.(f)))
            local_msg([ level, '.' f ' does not have the same class']);
        else
            switch class(a.(f))
                case 'struct'
                    local_addmsg(level, ['.' f ]);
                    local_compare(a.(f),b.(f),level)
                otherwise
                    aData = a.(f);
                    bData = b.(f);
                    isSizeOK = isequal(size(a.(f)),size(b.(f)));
                    isNumelOK = isequal(numel(a.(f)),numel(b.(f)));
                    isContentOK = false;
                    if isNumelOK
                        isContentOK = isequaln(a.(f)(:),b.(f)(:));
                    end
                    if isSizeOK && isContentOK
                        % do nothing - all good
                        local_addmsg(level, ['.' f '  - correct']);
                    elseif ~isContentOK
                        local_addmsg(level, ['.' f ' - &&&ERROR&&& has different CONTENT']);
                    else
                        local_addmsg(level, ['.' f ' - &&&ERROR&&& has different DIMENSIONS but similar content']);
                    end
            end
        end
    end
end

function msgAll = local_addmsg(level, msg)
    persistent msgStore
    
    verbose = true;
    switch msg
        case 'reset'
            msgStore = {};
        case 'retrieve'
            msgAll = msgStore;
        otherwise
            indent = repmat('  ',1,level);
            if verbose
                disp([indent,msg])
            end
            msgStore = [msgStore; [indent,msg]];
    end
    msgAll = msgStore;
end
