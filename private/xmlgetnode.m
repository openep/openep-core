function result = xmlgetnode(parentNode, searchName, level, searchType)
% XMLGETNODE finds a node with the specified name.
% Usage:
%   result = xmlgetnode(parentNode, searchName, level, searchType)
% Where:
%   parentNode is the parent xml node, or an xml filename
%   searchName is the name of the node to search for
%   level is the number of levels into the tree to look (pu Inf for all)
%   searchType can be 'first' or 'all'

% Author: Nick Linton (2013)
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------

% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------


if ischar(parentNode)
    %assume it is a filename
    parentNode = xmlread(parentNode);
end

switch searchType
    case 'first'
        findAll = false;
    case 'all'
        findAll = true;
end

result = sub_findnode(searchName, parentNode, level, findAll);

end

function result = sub_findnode(name, node, level, findAll)
    if strcmpi(name, char(node.getNodeName))
        result = node;
        return
    elseif level<=0
        result = [];
        return
    else
        cNodes = node.getChildNodes;
        nChild = cNodes.getLength;
        result = [];
        for i = 1:nChild
            newResult = sub_findnode(name, cNodes.item(i-1), level-1, findAll);
            if ~isempty(newResult)
                if findAll
                    result = [result , newResult]; %#ok<AGROW>
                else
                    result = newResult;
                    return
                end
            end
        end
    end
end
            