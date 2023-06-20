classdef LinearSolution < handle
    properties
        
        feas logical % boolean indicating feasibility
        bounded logical % boolean indicating bounded-ness
        x double % primal solution vector
        s double % dual solution vector
        undetermined logical % boolean vector indicating undetermined variables
        primal_zero logical % boolean vector indicating x_i = 0
        dual_zero logical % boolean vector indicating s_i = 0
        value double % optimal value
        
    end
    
    methods
        
        function obj = LinearSolution(feasible,bounded,x,s,value,undetermined,primal_zero,dual_zero)
            if ~feasible
                obj.feas = false;
            else
                obj.feas = true;
                if ~bounded
                    obj.bounded = false;
                else
                    obj.bounded = true;
                    obj.x = x;
                    obj.s = s;
                    obj.value = value;
                    if nargin > 5
                        obj.undetermined = undetermined;
                        obj.primal_zero = primal_zero;
                        obj.dual_zero = dual_zero;
                    end
                end
            end
        end
        
    end
    
end