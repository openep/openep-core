function [lowVArea, voltages, iTri, tr2] = getLowVoltageArea( userdata, varargin )
% GETLOWVOLTAGEAREA Returns the low voltage area
%
% Usage:
%   [lowVArea, voltages, iTri, tr2] = getLowVoltageArea( userdata, varargin )
% Where:
%   userdata  - see importcarto_mem
%   lowVArea  - the low voltage area (cm^2)
%   voltages - the voltages point values used to calculate areas
%   iTri - indexes into userdata.surface.triRep.Triangulation and refers
%          to the triangles that have voltage values within the range,
%          threshld
%   tr2 - a triangulation of all the triangles referenced in iTri.
%
%
% GETLOWVOLTAGEAREA accepts the following parameter-value pairs
%   'method'    {'map'} | 'egm'
%   'type'      {'bip'} | 'uni'
%   'threshold' {[0.0 0.5]} | array
%
% GETLOWVOLTAGEAREA Returns the surface area of the chamber with voltage
% less than the specified threshold, 0.5mV by default. By default low
% voltage area is calculated using the surface data (stored in
% userdat.surface). If 'method' is set to 'egm' then the bipolar voltage is 
% first interpolated from the bipolar electrogram data (stored in
% userdata.electric). If 'type' is set to 'uni' then unipolar voltages are
% used for surface area calculation.
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% [lowVoltageArea, voltages, iTri, tr2] = getLowVoltageArea(userdata, 'method', 'egm')
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 1; % UPDATE VALUE
method = 'map';
type = 'bip';
threshold = [0 0.5];
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'method'
                method = varargin{i+1};
            case 'type'
                type = varargin{i+1};
            case 'threshold'
                threshold = varargin{i+1};
        end
    end
end

switch lower(method)
    case 'map'
        switch lower(type)
            case 'bip'
                voltages = userdata.surface.act_bip(:,2);
            case 'uni'
                voltages = userdata.surface.uni_imp_frc(:,1);
        end
    case 'egm'
        switch lower(type)
            case 'bip'
                voltages = generateInterpData(userdata, 'bip-map');
            case 'uni'
                voltages = generateInterpData(userdata, 'uni-map');
        end
end

tr = userdata.surface.triRep;
[lowVArea, iTri, tr2] = local_calculateArea(tr, voltages, threshold);

    function [a, iTri, tr2] = local_calculateArea(tr, sI, threshold)
        sIFace = trVertToFaceData(tr, sI);
        iTri = zeros(size(sIFace));
        iTri(sIFace<threshold(2) & sIFace>threshold(1)) = 1; %keep these ones
        triangleInclude = tr.Triangulation;
        triangleInclude(~logical(iTri),:) = [];
        if ~isempty(triangleInclude)
            tr2 = TriRep(triangleInclude, tr.X);
            
            % Calculate the thresholded area
            areas2 = triarea(tr2);
            a = sum(areas2)/100;
        else
            a = 0;
            tr2 = [];
        end
    end

end

