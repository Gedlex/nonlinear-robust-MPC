% -------------------------------------------------------------------------
% File: PlanarQuadrotor.m
% Author: Alexander Erdin (aerdin@ethz.ch)
% Date: 21th April 2024
% License: MIT
% Reference:
%
% -------------------------------------------------------------------------
classdef PlanarQuadrotorPreStabilized < NonlinearSystem
    
    properties (GetAccess=public,SetAccess=protected)
        nx = 6              % Number of states
        nu = 2              % Number of inputs
        np = 1              % Number of parametric uncertainties
        nw = 1              % Number of disturbance states
        dt                  % Sampling interval

        K                   % Pre-stabilizing feedback gain
        
        g = 9.81            % [m/s^2] gravity
        l = 0.25            % [m] half-width of quad rotor
        m = 0.4811          % [kg] true mass of the quad rotor
        m_nom = 0.486       % [kg] nominal mass of the quad rotor
        J = 0.00383         % [kgm^2] moment of inertia
        
        plot = PlanarQuadrotorPlots() % Class for plotting
    end
    
    methods
        % state : x, z, phi, x_dot, z_dot, phi_dot
        function obj = initialize(obj,dt,options)
            arguments
                obj (1,1) PlanarQuadrotorPreStabilized
                dt (1,1) {mustBeNumeric, mustBePositive}
                options.K (:,:) {mustBeNumeric} = zeros(obj.nu,obj.nx)
                options.integrator (1,1) string {mustBeMember(options.integrator,["euler","euler_multi","rk4","rkdp"])} = 'rk4'
            end
            
            % Set sampling interval
            obj.dt = dt;

            % Pre-stabilizing feedback gain
            obj.K = options.K;
            
            % Define integrator
            obj.integrator = options.integrator;
        end
        
        function dt = f_fcn(obj,x) % State dynamics
            dt = [x(4,:).*obj.cosx(x(3,:)) - x(5,:).*obj.sinx(x(3,:)); % px
                  x(4,:).*obj.sinx(x(3,:)) + x(5,:).*obj.cosx(x(3,:)); % pz
                  x(6,:);                                              % phi
                  x(6,:).*x(5,:) - obj.g*obj.sinx(x(3,:));             % vx
                 -x(6,:).*x(4,:) - obj.g*obj.cosx(x(3,:));             % vz
                  zeros(1,size(x,2))];

            % Add pre-stabilizing feedback gain
            dt = dt + obj.B_fcn(x)*obj.K*x;
        end
        
        function dt = B_fcn(obj,x) % Input dynamics
            dt = [zeros(4,2);
                  1/obj.m_nom,  1/obj.m_nom;                           % vz
                  obj.l/obj.J, -obj.l/obj.J];                          % phi_dot
        end

        function dt = G_fcn(obj,x,u) % Parameter dynamics of the continuouse time system
            dt = [zeros(4,1);
                  u(1) + u(2);
                  0];
        end
        
        function dt = E_fcn(obj,x) % Disturbance dynamics
            dt = [0;
                  0;
                  0;
                  obj.cosx(x(3));
                 -obj.sinx(x(3));
                  0];
        end
    end
end