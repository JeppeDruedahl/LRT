classdef abb_par 
    
    properties 
        Param  = []; 
        b1     = []; 
        bL     = []; 
    end
    
    methods 
        function obj = abb_par(Param, b1, bL)
            if nargin>0
                obj.Param = Param; % parameters indexing the hermite basis functions
                obj.b1 = b1; 
                obj.bL = bL;
            end
        end

        function print(obj)
            fprintf('b1 = %8.4g\n', obj.b1); 
            fprintf('bL = %8.4g\n', obj.bL); 
            fprintf('Param(ibasis, jtau) = \n'); 
            disp(obj.Param)
        end
    end
end