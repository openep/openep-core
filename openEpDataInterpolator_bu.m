classdef openEpDataInterpolator_bu
    % OPENEPDATAINTERPOLATOR Creates objects for performing spatial
    % interpolation for OpenEP data
    %
    % Usage (Constructor):
    %   int = openEpDataInterpolator()
    %   int = openEpDataInterpolator(method)
    %   int = openEpDataInterpolator(method, options)
    %
    % Usage (interpolator):
    %   interpolated_function_at_query_points = int.interpolate(points,
    %       function_at_points, query_points);
    %
    % Where:
    %   method  - the method to use for interpolation, can be
    %             'scatteredInterpolant', 'localSmoothing' or 'radialBasis'
    %   options - a structure containing options pertaining to each
    %             interpolation method. Options which are not relevant to a
    %             particular interpolation method will be ignored.
    %       .interMethod
    %           - The interpolation method to use with
    %           scatteredInterpolant. See `help scatteredInterpolant`
    %       .exterMethod
    %           - The extrapolation method to use with
    %           scatteredInterpolant. See `help scatteredInterpolant`
    %       .smoothingLength
    %           - The smoothing length to use with localSmoothing.
    %           Consider a value in the range 5-10. Large values may overly
    %           smooth the resulting field, and small values may not
    %           provide enough coverage.
    %       .fillWith
    %           - The value to assign to the field at query points which
    %           fall outside the smoothingLength radius. Possible values
    %           are "nearest" or NaN. "nearest" performs nearest neighbour
    %           interpolation for query points outside the smoothingLength.
    %       .distanceThreshold
    %           - The distance threshold from known data to truncate the
    %           interpolated data.
    %
    % OPENEPDATAINTERPOLATOR Instances of openEpDataInterpolator perform
    % spatial interpolation. When instantiated, objects of this type
    % describe the interpolation scheme to be used. When ready to perform
    % interpolation, the `interpolation(...)` function should be called.
    %
    % Author: Steven Williams (2021) (Copyright)
    % SPDX-License-Identifier: Apache-2.0
    %
    % Modifications -
    %
    % Info on Code Testing:
    % ---------------------------------------------------------------
    % (1) Perform interpolation of electrogram voltage data using
    %     localSmoothing with default options.
    %
    % int = openEpDataInterpolator('localSmoothing');
    % egmX = getElectrogramX(userdata);
    % bip = getVoltages(userdata);
    % vtx = getVertices(userdata, 'used', false);
    % vertexVoltageData = int.interpolate(egmX, bip, vtx);
    %
    % (2) Visualise the above data
    %
    % figure;histogram(vertexVoltageData);
    % figure;drawMap(userdata, 'type', 'bip', 'coloraxis', [0.05 2])
    % title('Carto voltage data')
    % figure;drawMap(userdata, 'type', 'bip', 'coloraxis', [0.05 2], 'data', vertexVoltageData)
    % title('OpenEP voltage data')
    % ---------------------------------------------------------------
    %
    % ---------------------------------------------------------------
    % code
    % ---------------------------------------------------------------
    
    % Properties
    properties
        method
        interMethod
        exterMethod
        distanceThreshold
        smoothingLength
        fillWith
        rbfConstant
    end
    
    methods
        % Contructor
        function int = openEpDataInterpolator(varargin)
            
            % specify default values
            int.method = 'scatteredinterpolant';
            int.interMethod = 'linear';
            int.exterMethod = 'nearest';
            int.distanceThreshold = Inf;
            int.smoothingLength = 10;
            int.fillWith = 'nearest';
            int.rbfConstant = 1;
            
            % parse input data and warn if conflicts
            int.method = varargin{1};
            if nargin==2
                options = varargin{2};
                if isfield(options, 'interMethod')
                    if ~isempty(options.interMethod)
                        if ~strcmpi(int.method, 'scatteredInterpolant')
                            warning(['OPENEP/OPENEPDATAINTERPOLATOR: Setting interMethod property has no effect for interpolation method ' int.method]);
                        end
                        int.interMethod = options.interMethod;
                    end
                end
                if isfield(options, 'exterMethod')
                    if ~isempty(options.exterMethod)
                        if ~strcmpi(int.method, 'scatteredInterpolant')
                            warning(['OPENEP/OPENEPDATAINTERPOLATOR: Setting exterMethod property has no effect for interpolation method ' int.method]);
                        end
                        int.exterMethod = options.exterMethod;
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
                        end
                        int.smoothingLength = options.smoothingLength;
                    end
                end
                if isfield(options, 'fillWith')
                    if ~isempty(options.fillWith)
                        if ~strcmpi(int.method, 'localSmoothing')
                            warning(['OPENEP/OPENEPDATAINTERPOLATOR: Setting fillWith property has no effect for interpolation method ' int.method])
                        end
                        int.fillWith = options.fillWith;
                    end
                end
                if isfield(options, 'rbfConstant')
                    if ~isempty(options.fillWith)
                        if ~strcmpi(int.method, 'radialBasis')
                            warning(['OPENEP/OPENEPDATAINTERPOLATOR: Setting rbfConstant property has no effect for interpolation method ' int.method])
                        end
                        int.rbfConstant = options.rbfConstant;
                    end
                end
            end
        end
        
        % Perform Interpolation
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
                    
                    if isempty(x0) || isempty(x1) || isempty(d0)
                        error('OPENEP/OPENEPDATAINTERPOLATOR: Missing query or data points!');
                    end
                    % interpolation
                    [d1, ~] = localSmoothing(x0, d0, x1, int.smoothingLength, int.fillWith);
                    
                case 'radialbasis'
                    
                    if isempty(x0) || isempty(x1) || isempty(d0)
                        error('OPENEP/OPENEPDATAINTERPOLATOR: Missing query or data points!');
                    end
                    
                    % interpolation
                    op = rbfcreate(x0', d0', 'RBFFunction', 'multiquadric', 'RBFConstant', int.rbfConstant);
                    rbfcheck(op); %if this returns anything >10^-9 then consider increasing 'RBFConstant'
                    [d1, ~] = rbfinterp(x1', op);
                    d1 = d1';
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