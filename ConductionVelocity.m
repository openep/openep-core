classdef ConductionVelocity

    properties
        method = CVTriangulation;
        parameters;
    end

    methods

        function obj = ConductionVelocity(method)

           try
               obj.method = method;
           catch
           end
        
        end


        function run(obj, userdata)
            obj.method.run(userdata)
        end

    end

end
