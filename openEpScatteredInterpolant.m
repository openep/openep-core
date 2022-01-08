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
        function f_q = interpolate(obj, x, f_x, q)
            
            F = scatteredInterpolant(x(:,1), x(:,2), x(:,3), f_x, obj.interMethod, obj.exterMethod);
            f_q = F(q);
            
        end
    end
    
    % Getter and setter methods
    methods

       function obj = set.interMethod(obj, value)
          if ismember(value, {'linear', 'nearest', 'natural'})
              obj.interMethod = value;
          else
              error('Interpolation method must be one of ''linear'', ''nearest'' or ''natural''');
          end
       end

       function value = get.interMethod(obj)
          value = obj.interMethod;
       end
       
       function obj = set.exterMethod(obj, value)
           if ismember(value, {'linear', 'nearest', 'none'})
              obj.exterMethod = value;
          else
              error('Extrapolation method must be one of ''linear'', ''nearest'' or or ''none''');
          end
       end
       
       function value = get.exterMethod(obj)
          value = obj.exterMethod;
       end

    end
    
end