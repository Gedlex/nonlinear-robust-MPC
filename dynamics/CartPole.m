% -------------------------------------------------------------------------
% File: CaptureStabilization.m
% Author: Alexander Erdin (aerdin@ethz.ch)
% Date: 21th April 2024
% License: MIT
% Reference:
%
% -------------------------------------------------------------------------
classdef CartPole < NonlinearSystem
    
    properties (GetAccess=public,SetAccess=protected)
        nx = 4              % Number of states
        nu = 1              % Number of inputs
        np = 1              % Number of parametric uncertainties
        nw = 2              % Number of disturbance states
        dt                  % Sampling interval
        
        g = 9.81            % [m/s^2] gravity
        l_nom = 1           % [kgm^2] nominal half-length of the pole
        m_p = 1             % [kg] mass of the pole
        m_c = 1             % [kg] mass of the cart
        
        plot = PlanarQuadrotorPlots(xLabels={'x','\theta','{\dot x}','{\dot \theta}'},uLabels={'F'},mapping={[1,3],[2,4]}) 
    end
    
    methods
        % state : x, theta, x_dot, theta_dot
        function obj = initialize(obj,dt,options)
            arguments
                obj (1,1) CartPole
                dt (1,1) {mustBeNumeric, mustBePositive}
                options.integrator (1,1) string {mustBeMember(options.integrator,["euler","euler_multi","rk4","rkdp"])} = 'multi'
            end

            % Set sampling interval
            obj.dt = dt;

            % Define integrator
            obj.integrator = options.integrator;
        end
        
        function dt = f_fcn(obj,x) % State dynamics
            m_tot = obj.m_c + obj.m_p;
            dd_theta = sin(x(2,:))*(m_tot*obj.g - obj.m_p*obj.l_nom*x(4,:)^2*cos(x(2,:))) / ...
                                   (obj.l_nom*(4/3*m_tot     -    obj.m_p*cos(x(2,:))^2));
            dd_x = obj.m_p*obj.l_nom*(x(4,:)^2*sin(x(2,:)) - dd_theta*cos(x(2,:))) / m_tot;
            dt = [x(3,:);
                  x(4,:);
                  dd_x;
                  dd_theta];
        end
        
        function dt = B_fcn(obj,x) % Input dynamics
            m_tot = obj.m_c + obj.m_p;
            dt = [0;
                  0;
                  1/m_tot;
                 -cos(x(2,:)) / (obj.l_nom*(4/3*m_tot - obj.m_p*cos(x(2,:))^2))];
        end
        
        function dt = G_fcn(obj,x,u) % Parametric uncertainty dynamics        
            dt = zeros(4,1);
        end
        
        function dt = E_fcn(obj,x) % Disturbance dynamics
            dt = [zeros(2);
                 eye(2)];
        end
    end
end