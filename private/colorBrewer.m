function rgb = colorBrewer(name)
% COLORBREWER Returns RGB data for the nice colors
%
% Usage:
%   rgb = colorBrewer(name)
% Where:
%   name - is 'r'|'g'|'b'|'p'|'o'|'y'
%   rgb - is the RGB colorspec
%
% COLORBREWER See http://colorbrewer2.org/
%
% Author: Steven Williams (2014)
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

if isnumeric(name)
    if name==1
        rgb = [228 26 28]/256;
    elseif name==2
        rgb = [55 126 184]/256;
    elseif name==3
        rgb = [77 175 74]/256;
    elseif name==4
        rgb = [152 78 163]/256;
    elseif name==5
        rgb = [255 127 0]/256;
    elseif name==6
        rgb = [255 255 51]/256;
    elseif name==7
        rgb = [166 86 40]/256;
    elseif name==8
        rgb = [247 129 191]/256;
    elseif name==9
        rgb = [153 153 153]/256;
    end
else
    switch name
        case 'r'
            rgb = [228 26 28]/256;
        case 'g'
            rgb = [77 175 74]/256;
        case 'b'
            rgb = [55 126 184]/256;
        case 'p'
            rgb = [152 78 163]/256;
        case 'o'
            rgb = [255 127 0]/256;
        case 'y'
            rgb = [255 255 51]/256;
        case 'g1'
            rgb = [27 158 119]/256;
        case 'o1'
            rgb = [217 95 2]/256;
        case 'p1'
            rgb = [117 112 179]/256;
    end
end
