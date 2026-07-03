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
%   'mapped'     {'all'}|'onlymapped'
%
% GETAREA Returns the surface area of an anatomical model. The anatomical
% model can first be closed (filling any holes) by specifying the 'method',
% 'fill' ('nofill' by default). Using the 'mapped' parameter specify
% whether to return the total area or only the area which has been mapped.
% The area which has been mapped is identified by having non-NaN values in
% userdata.surface.act_bip(:,2).
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
mapped = 'all';
if nargin > nStandardArgs
    for i = 1:2:nargin-1
        switch varargin{i}
            case 'method'
                method = varargin{i+1};
            case 'mapped'
                mapped = varargin{i+1};
        end
    end
end

switch method
    case "nofill"
        tr = getMesh(userdata);
        area = sum(real(triarea(tr))/100);
    case "fill"
        tr = getClosedSurface(userdata);
        area = sum(real(triarea(tr))/100); 
end

switch mapped
    case 'onlymapped'
        % subtract the area of the subset of triangles which has not been mapped
        sIFace = trVertToFaceData(tr, userdata.surface.act_bip(:,2));
        iTri = zeros(size(sIFace));
        iTri(isnan(sIFace)) = 1;
        triangleInclude = tr.Triangulation;
        triangleInclude(~logical(iTri),:) = [];
        if ~isempty(triangleInclude)
            tr2 = TriRep(double(triangleInclude), tr.X);
            areas2 = sum(real(triarea(tr2))/100);
        else
            areas2 = 0;
        end
        area = area - areas2;
    case 'all'
        % make no modifications
end

end