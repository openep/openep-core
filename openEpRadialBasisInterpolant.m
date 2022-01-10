classdef openEpRadialBasisInterpolant
    % OPENEPRADIALBASISINTERPOLANT Perform radial basis function
    % interpolation
    %
    % Usage:
    %   int = openEpDataInterpolator('radialbasis')
    %   [f_q, df_q] = int.interpolate(x, f_x, q);
    %
    % Where:
    %   x       - location points of data
    %   data    - data
    %   q       - location of query points
    %
    % OPENEPRADIALBASISINTERPOLANT is used for performing radial basis
    % function interpolation. Instantiate an object of the
    % OPENEEPDATAINTERPOLATOR class, setting `'method'` to 'radialbasis'.
    % The following properties may then be set:
    %       .basisFunction
    %           - TODO
    %       .doOptimisation
    %           - TODO
    %       .shapeParameter
    %           - TODO
    %       .distanceThreshold
    %           - The distance threshold from known data to truncate the
    %           interpolated data.
    %
    % Author: Steven Williams / Chris O'Shea / Adam Connolly(2022) (Copyright)
    % SPDX-License-Identifier: Apache-2.0
    %
    % Modifications -
    %
    % See also: openEpDataInterpolator, openEpLocalSmoothing, openEpScatteredInterpolant
    %
    % Info on Code Testing:
    % ---------------------------------------------------------------
    % 
    % ---------------------------------------------------------------
    %
    % ---------------------------------------------------------------
    % code
    % ---------------------------------------------------------------

    properties
        basisFunction = 'multiquadric';
        doOptimisation = 'false';
        shapeParameter = 1;
        distanceThreshold = 5;
    end

    methods

        % Main method
        function [f_q, df_q] = interpolate(obj, x, f_x, q)

            rbfoptions.basisFunction=obj.basisFunction;
            rbfoptions.shapeParameter=obj.shapeParameter;
            rbfoptions.doOptimisation=obj.doOptimisation;

            % remove any NaN values from data
            iNan = isnan(f_x);
            x(iNan,:) = [];
            f_x(iNan) = [];
            
            [f_q, df_q] = rbf_interpolator(f_x,x,q,rbfoptions);
            f_q=f_q';
            
            % filter by distance
            f_q = filterByDistance(f_q, q, x, obj.distanceThreshold);
            df_q = filterByDistance(df_q, q, x, obj.distanceThreshold);
        end
    end

    % Getter and setter methods
    methods
        function obj = set(obj, property, value)
            switch lower(property)
                case 'basisfunction'
                    if ismember(value, {'cubic', 'linear', 'thinplate', 'gaussian', 'multiquadric'})
                        obj.basisFunction = value;
                    else
                        error('Basis must be one of ''cubic'', ''linear'', ''thinplate'', ''gaussian'' or ''multiquadric''');
                    end
                case 'dooptimisation'
                    obj.doOptimisation = value;
                case 'shapeparameter'
                    obj.shapeParameter = value;
                case 'distancethreshold'
                    obj.distanceThreshold = value;
            end
        end

        function value=get.doOptimisation(obj)
            value=obj.doOptimisation;
        end

        function value=get.shapeParameter(obj)
            value=obj.shapeParameter;
        end
    end
end