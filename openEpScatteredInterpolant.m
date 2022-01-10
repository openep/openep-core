classdef openEpScatteredInterpolant < matlab.mixin.SetGet

    properties
        interMethod = 'linear';
        exterMethod = 'nearest';
        distanceThreshold = 10;
    end

    methods
        % Default implicit constructor.
        function obj = openEpScatteredInterpolantd(distanceThreshold)
            obj.distanceThreshold = distanceThreshold;
        end

        % Main method
        function [f_q, df_q] = interpolate(obj, x, f_x, q)

            % calculate the interpolated values
            F = scatteredInterpolant(x(:,1), x(:,2), x(:,3), f_x, obj.interMethod, obj.exterMethod);
            f_q = F(q);
            f_q = filterByDistance(f_q, q, x, obj.distanceThreshold);

            % calculate the gradient of these values
            shrinkFactor = 1;
            tri = boundary(q, shrinkFactor);
            tr = TriRep(tri,q(:,1),q(:,2),q(:,3));
            [~, df_q] = trigrad(tr, f_q);
            df_q = filterByDistance(df_q, q, x, obj.distanceThreshold);
        end
    end

    % Getter and setter methods
    methods

        function set.interMethod(obj, value)
            if ismember(value, {'linear', 'nearest', 'natural'})
                obj.interMethod = value;
            else
                error('Interpolation method must be one of ''linear'', ''nearest'' or ''natural''');
            end
        end

        function value = get.interMethod(obj)
            value = obj.interMethod;
        end

        function set.exterMethod(obj, value)
            if ismember(value, {'linear', 'nearest', 'none'})
                obj.exterMethod = value;
            else
                error('Extrapolation method must be one of ''linear'', ''nearest'' or or ''none''');
            end
        end

        function value = get.exterMethod(obj)
            value = obj.exterMethod;
        end

        function set.distanceThreshold(obj, value)
            if isnumeric(value)
                obj.distanceThreshold = value
            else
                error('OPENEP/OPENEPSCATTEREDINTERPOLANT: distanceThreshold must be numeric.')
            end
        end

        function value = get.distanceThreshold(obj)
            value = obj.distanceThreshold;
        end

    end

end