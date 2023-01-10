function area = getArea( userdata, varargin )
% GETAREA Returns the surface area of an anatomical model
%
% Usage:
%   area = getArea( userdata )
% Where:
%   userdata  - see importcarto_mem
%   area  - the surface area (cm^2)
%
% GETAREA accepts the following parameter-value pairs
%   'method'     {'nofill'}|'fill'
%
% GETAREA Returns the surface area of an anatomical model. The anatomical
% model can first be closed (filling any holes) by specifying the 'method',
% 'fill' ('nofill' by default).
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% area = getArea( userdata, 'method', 'fill' )
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 1; % UPDATE VALUE
method = 'nofill';
if nargin > nStandardArgs
    for i = 1:2:nargin-1
        switch varargin{i}
            case 'method'
                method = varargin{i+1};
        end
    end
end

switch method
    case "nofill"
        area = sum(real(triarea(getMesh(userdata)))/100);
    case "fill"
        tr = getClosedSurface(userdata);
        area = sum(real(triarea(tr))/100); 
end
end
