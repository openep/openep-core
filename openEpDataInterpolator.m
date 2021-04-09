classdef openEpDataInterpolator
    % OPENEPDATAINTERPOLATOR Performs data interpolation for OpenEP
    
    properties
        method
    end
    
    properties (Access = private)
        interMethod
        exterMethod
        distanceThreshold
        smoothingLength
        fillWith
    end
    
    methods
        % Contructor
        function int = openEpDataInterpolator(varargin)
            
            % specify default values
            int.interMethod = 'linear';
            int.exterMethod = 'linear';
            int.distanceThreshold = Inf;
            int.smoothingLength = 10;
            int.fillWith = NaN;
            
            int.method = 'scatteredinterpolant';
            
            if nargin==1
                int.method = varargin{1};
            end
            if nargin==2
                options = varargin{2};
                if isfield(options, 'interMethod')
                    if ~isempty(options.interMethod)
                        if ~strcmpi(int.method, 'scatteredInterpolant')
                            warning(['OPENEP/OPENEPDATAINTERPOLATOR: Setting interMethod property has no effect for interpolation method ' int.method]);
                            int.interMethod = options.interMethod;
                        end
                    end
                end
                if isfield(options, 'exterMethod')
                    if ~isempty(options.exterMethod)
                        if ~strcmpi(int.method, 'scatteredInterpolant')
                            warning(['OPENEP/OPENEPDATAINTERPOLATOR: Setting exterMethod property has no effect for interpolation method ' int.method]);
                            int.exterMethod = options.exterMethod;
                        end
                    end
                end
                if isfield(options, 'distanceThreshold')
                    if ~isempty(options.distanceThreshold)
                        int.distanceThreshold = options.distanceThreshold;
                    end
                end
                if isfield(options, 'smoothingLength')
                    if ~isempty(options.smoothingLength)
                        if ~strcmpi(int.method, 'localSmoothing')
                            warning(['OPENEP/OPENEPDATAINTERPOLATOR: Setting fillWith property has no effect for interpolation method ' int.method])
                            int.smoothingLength = options.smoothingLength;
                        end
                    end
                end
                if isfield(options, 'fillWith')
                    if ~isempty(options.fillWith)
                        if ~strcmpi(int.method, 'localSmoothing')
                            warning(['OPENEP/OPENEPDATAINTERPOLATOR: Setting fillWith property has no effect for interpolation method ' int.method])
                            int.fillWith = options.fillWith;
                        end
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
                    
                    % error checking
                    if numel(d0) < 4
                        error('OPENEP/OPENEPDATAINTERPOLATOR: Cannot perform interpolation with less than four data points');
                    end
                    
                    % interpolation
                    F = scatteredInterpolant(x0(:,1), x0(:,2), x0(:,3), d0, int.interMethod, int.exterMethod);
                    d1 = F(x1);
                    
                case 'localsmoothing'
                    
                    % error checking
                         % % % AC: is there any appropriate error checking to do
                    % here?
                    
                    % interpolation
                    d1 = localSmoothing(x0, d0, x1, int.smoothingLength, int.fillWith);
                    
                case 'radialbasis'
                    
                    % error checking
                    
                    % interpolation
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