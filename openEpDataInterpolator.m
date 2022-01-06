classdef openEpDataInterpolator
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
        distanceThreshold
    end
    
    properties (Access = private)
        call_interpolator
    end
    
    methods
        % Contructor
        function obj = openEpDataInterpolator(method)
                        
            possibleMethods = {'scatteredinterpolant' ...
                , 'localsmoothing' ...
                , 'radialbasis' ...
                };
            
            defaultMethod = 'scatteredinterpolant';
            
            obj.distanceThreshold = 10;
            
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
                                obj.call_interpolator = openEpLocalSmoothingInterpolant;
                            case 'radialbasis'
                                obj.call_interpolator = openEpRadialBasisInterpolant;
                        end
                    end
                end      
            end
                
        end
        
        function obj = set(obj, property, value)
          switch property
          case 'basisFunction'    
          if ismember(value, {'cubic', 'linear', 'thinplate', 'gaussian', 'multiquadric'})
              obj.basisFunction = value;
          else
              error('Basis must be one of ''cubic'', ''linear'', ''thinplate'', ''gaussian'' or ''multiquadric''');
          end
          case 'doOptimisation' 
              obj.doOptimisation = value;
          case 'shapeParameter'
              obj.shapeParameter = value;  
         end  
       end
        
        function f_q = interpolate(obj, x, f_x, q)
            
            f_q = obj.call_interpolator.interpolate(x, f_x, q);
            
        end
        
        function disp(obj)
            display(obj.call_interpolator);
        end
    end
end