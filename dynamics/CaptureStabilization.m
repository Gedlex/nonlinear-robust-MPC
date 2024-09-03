% -------------------------------------------------------------------------
% File: CaptureStabilization.m
% Author: Alexander Erdin (aerdin@ethz.ch)
% Date: 21th April 2024
% License: MIT
% Reference:
%
% -------------------------------------------------------------------------
classdef CaptureStabilization < NonlinearSystem
    
    properties (GetAccess=public,SetAccess=protected)
        nx = 6              % Number of states
        nu = 2              % Number of inputs
        np = 1              % Number of parametric uncertainties
        nw = 3              % Number of disturbance states
        dt                  % Sampling interval
        
        l = 1               % [m] half-width of quad rotor
        j_nom = 1           % [kgm^2] nominal inertia of the satellite
        m = 1               % [kg] mass of the satalite
        
        plot = PlanarQuadrotorPlots() % Class for plotting
    end
    
    methods
        % state : x, y, theta x_dot, y_dot, theta_dot
        function obj = initialize(obj,dt,options)
            arguments
                obj (1,1) CaptureStabilization
                dt (1,1) {mustBeNumeric, mustBePositive}
                options.integrator (1,1) string {mustBeMember(options.integrator,["euler","euler_multi","rk4","rkdp"])} = 'multi'
            end

            % Set sampling interval
            obj.dt = dt;

            % Define integrator
            obj.integrator = options.integrator;
        end
        
        function dt = f_fcn(obj,x) % State dynamics
            dt = [x(4,:);
                  x(5,:);
                  x(6,:);
                  0;
                  0;
                  0];
        end
        
        function dt = B_fcn(obj,x) % Input dynamics
            theta = x(3);
            dt = [0,                      0;
                  0,                      0;
                  0,                      0;
                  obj.cosx(theta)/obj.m, -obj.sinx(theta)/obj.m;
                  obj.sinx(theta)/obj.m,  obj.cosx(theta)/obj.m;
                  obj.l/obj.j_nom,        0];
        end
        
        function dt = G_fcn(obj,x,u) % Parametric uncertainty dynamics        
            dt = [zeros(5,1);
                  obj.l*u(1)];
        end
        
        function dt = E_fcn(obj,x) % Disturbance dynamics
            E = zeros(obj.nx,obj.nw);
            E(5,1)=1; E(4,2)=1; E(6,3)=1;
            dt = E;
        end
    end
end