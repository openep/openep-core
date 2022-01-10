classdef openEpDataInterpolator < matlab.mixin.SetGet
    % OPENEPDATAINTERPOLATOR Creates objects for performing spatial
    % interpolation for OpenEP data
    %
    % Usage:
    %   int = openEpDataInterpolator()
    %   int = openEpDataInterpolator(method)
    %   [f_q, df_q] = int.interpolate(x, f_x, q);
    %
    % Where:
    %   method  - the method to use for interpolation, can be
    %             'scatteredInterpolant', 'localSmoothing' or
    %             'radialBasis'
    %   x       - location points of data
    %   data    - data
    %   q       - location of query points
    %
    % OPENEPDATAINTERPOLATOR Instances of openEpDataInterpolator perform
    % spatial interpolation. When instantiated, objects of this type
    % describe the interpolation scheme to be used. When ready to perform
    % interpolation, the `interpolate(...)` function should be called.
    % Interpolation-scheme specific settings can be configured. See help
    % <scheme name> for details:
    %       help OPENEPRADIALBASISINTERPOLANT
    %       help OPENEPSCATTEREDINTERPOLANT
    %       help OPENEPLOCALSMOOTHING
    %
    % Author: Steven Williams / Adam Connolly (2022) (Copyright)
    % SPDX-License-Identifier: Apache-2.0
    %
    % Modifications -
    %
    % See also: openEpRadialBasisInterpolant, openEpLocalSmoothing, openEpScatteredInterpolant
    %
    % Info on Code Testing:
    % ---------------------------------------------------------------
    % (1) Perform interpolation of electrogram voltage data using default options.
    %
    % int = openEpDataInterpolator('localSmoothing');
    % egmX = getElectrogramX(userdata);
    % bip = getVoltages(userdata);
    % vtx = getVertices(userdata, 'used', false);
    % vertexVoltageData = int.interpolate(egmX, bip, vtx);
    %
    % (2) Visualise
    %
    % figure;histogram(vertexVoltageData);
    % figure;drawMap(userdata, 'type', 'bip', 'coloraxis', [0.05 2])
    % title('Original voltage data')
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
    end

    properties (Access = private)
        call_interpolator
    end

    properties (Dependent)
        interMethod
        exterMethod
        basisFunction
        doOptimisation
        shapeParameter
        smoothingLength
        fillWith
        distanceThreshold
    end

    methods
        % Contructor
        function obj = openEpDataInterpolator(method)

            possibleMethods = {'scatteredinterpolant' ...
                , 'localsmoothing' ...
                , 'radialbasis' ...
                };

            defaultMethod = 'scatteredinterpolant';

            if nargin == 0
                obj.method = defaultMethod;
                obj.call_interpolator = openEpScatteredInterpolant;
            else
                if nargin > 1
                    error('Too many input arguments!');
                else
                    if ~ismember(lower(method), possibleMethods)
                        error('Choose from ''scatteredinterpolant'', ''localsmoothing'' or ''radialbasis''');
                    else
                        obj.method = lower(method);
                        switch obj.method
                            case 'scatteredinterpolant'
                                obj.call_interpolator = openEpScatteredInterpolant;
                            case 'localsmoothing'
                                obj.call_interpolator = openEpLocalSmoothing;
                            case 'radialbasis'
                                obj.call_interpolator = openEpRadialBasisInterpolant;
                        end
                    end
                end
            end
        end

        % overloaded set method
        function obj = set(obj, property, value)
            set(obj.call_interpolator, property, value)
        end

        % overloaded get methods
        function value = get(obj, property)
            value = get(obj.call_interpolator, property);
        end

        % overloaded dot syntax
        function set.interMethod(obj, value)
            obj.call_interpolator.interMethod = value;
        end
        function set.exterMethod(obj, value)
            obj.call_interpolator.exterMethod = value;
        end
        function set.basisFunction(obj, value)
            obj.call_interpolator.basisFunction = value;
        end
        function set.doOptimisation(obj, value)
            obj.call_interpolator.doOptimisation = value;
        end
        function set.shapeParameter(obj, value)
            obj.call_interpolator.shapeParameter = value;
        end
        function set.smoothingLength(obj, value)
            obj.call_interpolator.smoothingLength = value;
        end
        function set.fillWith(obj, value)
            obj.call_interpolator.fillWith = value;
        end
        function set.distanceThreshold(obj, value)
            obj.call_interpolator.distanceThreshold = value;
        end

        % overloaded display method
        function disp(obj)
            display(obj.call_interpolator);
        end

        % main interpolation function
        function [f_q, df_q] = interpolate(obj, x, f_x, q)

            % preprocess to f_x and x to handle duplicate data points
            [x, f_x] = removeDuplicateDataPoints(x, f_x);

            % perform interpolation using the specified interpolation method
            [f_q, df_q] = obj.call_interpolator.interpolate(x, f_x, q);
        end
    end
end
