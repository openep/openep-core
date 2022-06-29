classdef ConductionVelocity < matlab.mixin.SetGet & matlab.System
    % ConductionVelocity Creates objects for calculating conduction
    % velocity from OpenEP datasets.
    %
    % Usage:
    %   cv = CondictionVelocity()
    %   cv = CondictionVelocity(method);
    %   [cv_values, cv_X] = cv.run(userdata);
    %
    % Where:
    %   method  - the method to use for calculating conduction velocity.
    %             Must be a class that implements a `run` method for
    %             calculating conduction velocity. Classes available
    %             in OpenEP:
    %               - CVTriangulation
    %               - CVMappingCosineFit
    %               - CVPlanarFit
    %   cv_values - calculated conduction velocity values
    %   cv_X      - location (3D coordinates) of each conduction velocity value.
    %
    % Instances of ConductionVelocity class calculate the conduction velocity
    % using a specific algorithm. When instantiated, objects of this type
    % describe the algorithm and parameters to be used. When ready to perform
    % the calculation, the `run(userdata)` function should be called.
    % Algorithm-specific parameters can be configured directly from the
    % ConductionVelocity instance. See help <class name> for details:
    %       help CVTriangulation
    %
    % Author: Paul Smith (2022) (Copyright)
    % SPDX-License-Identifier: Apache-2.0
    %
    % Info on Code Testing:
    % ---------------------------------------------------------------
    % (1) Perform conduction velocity calculation using the default
    % algorithm (CVTriangulation) and parameters.
    %
    % cv = ConductionVelocity();
    % [cv_values, cv_X] = cv.run(userdata);
    %
    % ---------------------------------------------------------------
    %
    % ---------------------------------------------------------------
    % code
    % ---------------------------------------------------------------

    properties
        method = CVTriangulation;
    end

    properties(Dependent = true)
        latDif;                     % CVTriangulation
        elecDis;                    % CVTriangulation
        maxCVallowed;               % CVPlanarFit
        minCVallowed;               % CVTriangulation
        minTheta;                   % CVTriangulation
        bipolarVoltageThreshold;    % CVPlanarFit
        plot;                       % CVTriangulation, CVMappingCosineFit
    end

    methods

        % Constructor
        function obj = ConductionVelocity(method)

            % If no method is passed, the default (CVTriangulation) will be
            % kept
            try
                obj.method = method;
            catch
            end

        end

        % overloaded set method
        function obj = set(obj, property, value)

            % TODO: correctly display set(obj)
            % See the description of `s = set(H)` at:
            % https://uk.mathworks.com/help/matlab/ref/set.html
            % Currently, if `set(obj)` is called on the ConductionVelocity class,
            % property and value will not exist and an error will be raised

            if strcmp(property, 'method')
                obj.method = value;
            else
                % Ensure ConductionVelocity and method class have same
                % properties
                set(obj.method, property, value);
            end
            
        end

        % overloaded get methods
        function value = get(obj, property)

            % TODO: correctly display get(obj)
            % See the description of `get(sys)` at:
            % https://uk.mathworks.com/help/control/ref/lti.get.html
            % Currently, if `get(obj)` is called on the ConductionVelocity class,
            % property will not exist

            if strcmp(property, 'method')
                value = obj.method;
            else
                value = get(obj.method, property);
            end

        end

        % overloaded dot syntax
        function set.latDif(obj, value)
            obj.method.latDif = value;
        end
        function value = get.latDif(obj)
            value = obj.method.latDif;
        end

        function set.elecDis(obj, value)
            obj.method.elecDis = value;
        end
        function value = get.elecDis(obj)
            value = obj.method.elecDis;
        end

        function set.maxCVallowed(obj, value)
            obj.method.maxCVallowed = value;
        end
        function value = get.maxCVallowed(obj)
            value = obj.method.maxCVallowed;
        end

        function set.minCVallowed(obj, value)
            obj.method.minCVallowed = value;
        end
        function value = get.minCVallowed(obj)
            value = obj.method.minCVallowed;
        end

        function set.minTheta(obj, value)
            obj.method.minTheta = value;
        end
        function value = get.minTheta(obj)
            value = obj.method.minTheta;
        end

        function set.bipolarVoltageThreshold(obj, value)
            obj.method.bipolarVoltageThreshold = value;
        end
        function value = get.bipolarVoltageThreshold(obj)
            value = obj.method.bipolarVoltageThreshold;
        end

        function set.plot(obj, value)
            obj.method.plot = logical(value);
        end
        function value = get.plot(obj)
            value = obj.method.plot;
        end

        % Run method
        function [cv_values, cv_X] = run(obj, userdata)
            [cv_values, cv_X] = obj.method.run(userdata);
        end

    end

    methods (Access = protected)

        % Determine which properties are active (for displaying info)
        function inactive = isInactivePropertyImpl(obj,propertyName)

            propNames = properties(obj.method);
            if strcmp('method', propertyName) || any(strcmp(propNames, propertyName))
                inactive = false;
            else
                inactive = true;
            end

        end

    end

end
