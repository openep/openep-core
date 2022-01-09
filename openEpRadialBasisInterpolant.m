classdef openEpRadialBasisInterpolant

    properties
        basisFunction = 'multiquadric';
        doOptimisation = 'false';
        shapeParameter = 1;
        distanceThreshold = 10;
    end

    methods

        % Main method
        function f_q = interpolate(obj, x, f_x, q)

            rbfoptions.basisFunction=obj.basisFunction;
            rbfoptions.shapeParameter=obj.shapeParameter;
            rbfoptions.doOptimisation=obj.doOptimisation;
            
            f_q = rbf_interpolator(f_x,x,q,rbfoptions);
            f_q=f_q';

            % accept only those interpolated values in proximity to actual values
            vtx = getVerticesNearMappingPoints(userdata, DISTANCETHRESHOLD);
            cv(~vtx) = [];
            cvX(~vtx,:) = [];
            n(~vtx,:) = [];
            u(~vtx,:) = [];
            disp(['OPENEP/GETCONDUCTIONVELOCITY: ' num2str(sum(~vtx)) ' CV values were removed which were more than ' num2str(DISTANCETHRESHOLD) 'mm from a mapping point']);

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