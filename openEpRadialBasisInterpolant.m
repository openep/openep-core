classdef openEpRadialBasisInterpolant
    
    properties
        basisFunction = 'multiquadric';
        doOptimisation = 'false';
        shapeParameter=1;
        distanceThreshold=10;
    end
    
    methods

        
        function f_q = interpolate(obj, x, f_x, q)
            rbfoptions.basisFunction=obj.basisFunction;
            rbfoptions.shapeParameter=obj.shapeParameter;
            rbfoptions.doOptimisation=obj.doOptimisation;
            f_q = rbf_interpolator(f_x,x,q,rbfoptions); 
            f_q=f_q';
        end
    end
    
  methods
       
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
          case 'distanceThreshold'
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