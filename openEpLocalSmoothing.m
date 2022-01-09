classdef openEpLocalSmoothing < matlab.mixin.SetGet
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

    properties
        smoothingLength = 5;
        fillWith = 'nearest';
        distanceThreshold = 10;
    end

    methods
        % Default implicit constructor.
        function obj = openEpScatteredInterpolantd(distanceThreshold)
            obj.distanceThreshold = distanceThreshold;
        end

        % Main method
        function [f_q, df_q] = interpolate(obj, x, f_x, q)

            [f_q, df_q] = localSmoothing(x, f_x, q, obj.smoothingLength, obj.fillWith);

            f_q = filterByDistance(f_q, q, x, obj.distanceThreshold);
            df_q = filterByDistance(df_q, q, x, obj.distanceThreshold);

        end
    end

    % Getter and setter methods
    methods

        function set.smoothingLength(obj, value)
            if isnumeric(value)
                obj.smoothingLength = value;
            else
                error('OPENEP/OPENEPLOCALSMOOTHING: The value for smoothingLength must be numeric.');
            end
        end

        function value = get.smoothingLength(obj)
            value = obj.smoothingLength;
        end

        function set.fillWith(obj, value)
            if isnan(value)
                obj.fillWith = value;
            elseif ischar && strcmpi(value, 'nearest')
                obj.fillWith = value;
            else
                error('OPENEP/OPENEPLOCALSMOOTHING: THe value for fillWith must be "nearest" or NaN.')
            end
        end

        function value = get.fillWith(obj)
            value = obj.fillWith;
        end

    end

end