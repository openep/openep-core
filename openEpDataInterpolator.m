classdef openEpDataInterpolator
    % OPENEPDATAINTERPOLATOR Performs data interpolation for OpenEP
    
    properties
        method
        interMethod
        exterMethod
        distanceThreshold
        smoothingLength
    end
    
    methods
        % Contructor
        function int = openEpDataInterpolator(varargin)
            
            % specify default values
            int.interMethod = 'linear';
            int.exterMethod = 'linear';
            int.distanceThreshold = Inf;
            int.smoothingLength = 10;
            int.method = 'scatteredinterpolant';
            
            if nargin==1
                int.method = varargin{1};
            end
            if nargin==2
                options = varargin{2};
                if isfield(options, 'interMethod')
                    if ~isempty(options.interMethod)
                        int.interMethod = options.interMethod;
                    end
                end
                if isfield(options, 'exterMethod')
                    if ~isempty(options.exterMethod)
                        int.exterMethod = options.exterMethod;
                    end
                end
                if isfield(options, 'distanceThreshold')
                    if ~isempty(options.distanceThreshold)
                        int.distanceThreshold = options.distanceThreshold;
                    end
                end
            end
        end
        
        % Interpolation
        function d1 = interpolate(int, x0,d0,x1)
            
            % remove any data with NaNs
            tempData = d0;
            tempCoords = x0;
            iNaN = isnan(tempData);
            tempData(iNaN) = [];
            tempCoords(iNaN,:) = [];
            d0 = tempData;
            x0 = tempCoords;
            
            % perform interpolation
            switch lower(int.method)
                case 'scatteredinterpolant'
                    if numel(d0) < 4
                        error('OPENEP/OPENEPDATAINTERPOLATOR: Cannot perform interpolation with less than four data points');
                    end
                    F = scatteredInterpolant(x0(:,1), x0(:,2), x0(:,3) ...
                        , d0, int.interMethod, int.exterMethod);
                    d1 = F(x1);
                    
                case 'localsmoothing'
                    
                case 'radialbasis'
            end
            
            % remove interpolated data that is too far from real data
            id = knnsearch(x0, x1);
            cPts = x0(id,:); %c for closest
            d = distBetweenPoints(cPts, x1);
            thresholdDistance = zeros(size(d));
            thresholdDistance(d>int.distanceThreshold) = 1;
            nanset = logical(thresholdDistance);
            d1(nanset) = NaN;
        end
    end
end