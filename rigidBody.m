classdef rigidBody
    % rigidBody class
    properties
        D       % dimension
        R       % radius for visualization
        M       % scalar or D x D Mass matrix
        cm      % center of mass
        J       % D x D moment of inertia matrix
        N       % number of ports
        x       % N x D array of ports positions
    end
    
    methods
        function obj = rigidBody(D,R,M,cm,J,x)
            % constructor of an object part of an embedded 1D fiber
            % in:       D = dimension 2 or 3
            %           R = object radius for visualization
            %           M = mass matrix
            %           c = center of mass
            %           J = moment of inertia matrix
            %           ports = N x D array of ports positions 
            % out: rigidBody object
            
            obj.D  = D;
            obj.R  = R;
            obj.N  = length(x);
            obj.M  = M;
            obj.cm = cm;
            obj.J  = J;
            obj.x = x;  
        end
    end
end