% -------------------------------------------------------------------------
% File: Parameters.m
% Author: Alexander Erdin (aerdin@ethz.ch)
% Date: 21th April 2024
% License: MIT
% Reference:
%
% -------------------------------------------------------------------------
classdef Parameters
    
    properties
        label (1,1) string                                                 % Label for the parameters
        version (1,1) {mustBeInteger,mustBePositive} = 1                   % Version of the parameters
        
        N  (1,1) {mustBeInteger,mustBeNonnegative}                         % Horizon (prediction steps)
        dt (1,1) {mustBeNumeric}                                           % Sampling interval
        
        x0 (:,:) {mustBeNumeric,mustBeVector(x0,"allow-all-empties")}      % Initial state
        x_ref (:,:) {mustBeNumeric,mustBeVector(x_ref,"allow-all-empties")}% Reference state
        u_ref (:,:) {mustBeA(u_ref,["function_handle","double"])}          % Reference input
        
        F_u (:,:) {mustBeNumeric}                                          % Input polytope
        b_u (:,:) {mustBeNumeric,mustBeVector(b_u,"allow-all-empties")}    % Input back-off
        F_x (:,:) {mustBeNumeric}                                          % State polytope
        b_x (:,:) {mustBeNumeric,mustBeVector(b_x,"allow-all-empties")}    % State back-off
        h_obs (1,1) {mustBeA(h_obs,["function_handle","double"])} = NaN    % Nonlinear obstacle constraint
        
        nw (1,1) = NaN                                                     % Number of disturbance states
        w_max (1,1) = NaN                                                  % Maximum disturbance
        theta_v (:,2) {mustBeNumeric} = NaN                                % Compact parameter set
        theta_true (:,:) {mustBeNumeric,mustBeVector(theta_true,"allow-all-empties")}
        
        Q_cost (:,:) {mustBeNumeric}                                       % State cost
        R_cost (:,:) {mustBeNumeric}                                       % Input cost
        P_cost (:,:) {mustBeNumeric}                                       % Terminal cost
        regularizer (1,1) {mustBeNonnegative}                              % Regularizer
        
        integrator (1,1) string {mustBeMember(integrator,["euler","euler_multi","rk4","rkdp"])} = "euler"
    end
    
    properties (Dependent)
        T  (1,1) {mustBeNumeric}                                           % Time horizon
        F_d (:,:) {mustBeNumeric}                                          % Disturbance polytope
        b_d (:,:) {mustBeNumeric,mustBeVector(b_d,"allow-all-empties")}    % Disturbance back-off
        F_theta (:,:) {mustBeNumeric}                                      % Theta polytope
        b_theta (:,:) {mustBeNumeric,mustBeVector(b_theta,"allow-all-empties")} % Theta back-off
        disturbance_v (:,:) {mustBeNumeric}                                % Compact disturbance set
        param_uncertainty (1,1) logical                                    % Boolean for parametric uncertainty
    end
    
    methods
        function obj = Parameters(params)
            switch nargin
                case 0
                    return
                case 1
                    % Get struct fieldnames
			        fields = fieldnames(params);

                    % Loop through struct and search for matching entry
                    for i = 1:numel(fields)
                        % Check if property exists
                        if ~isprop(obj, fields{i})
                            id = 'Parameters:InvalidProperty';
                            error(id,'%s is not a valid property',fields{i})
                        end
                        
                        % Set value
                        obj.(fields{i}) = params.(fields{i});
                    end
                otherwise
                    error('Too many input arguments.')
            end
        end
        
        function obj = set.F_x(obj,F)
            obj.checkPolytopes(F)
            obj.F_x = F;
        end
        
        function obj = set.F_u(obj,F)
            obj.checkPolytopes(F)
            obj.F_u = F;
        end
        
        function obj = set.nw(obj,n)
            mustBeNonnegative(n);
            mustBeInteger(n);
            obj.nw = n;
        end
        
        function obj = set.w_max(obj,w)
            mustBeNonnegative(w);
            mustBeNumeric(w)
            obj.w_max = w;
        end
        
        function obj = set.theta_v(obj,theta)
            if theta(:,1) > 0 ||  theta(:,2) < 0
                id = 'Parameters:EmptySet';
                error(id,['Theta must span a nonempty set containing the origin. ' ...
                          'The first column indicates the minimum and the second the maximum.'])
            end
            obj.theta_v = theta;
        end
        
        function obj = set.Q_cost(obj,Q)
            obj.mustBePositiveDefinite(Q)
            obj.Q_cost = Q;
        end
        
        function obj = set.R_cost(obj,R)
            obj.mustBePositiveDefinite(R)
            obj.R_cost = R;
        end
        
        function obj = set.P_cost(obj,P)
            obj.mustBePositiveDefinite(P)
            obj.P_cost = P;
        end
        
        function T = get.T(obj)
            T = obj.dt*obj.N;
        end
        
        function F = get.F_d(obj)
            obj.mustBeSet(obj.nw,'nw')
            F = [eye(obj.nw); -eye(obj.nw)];
        end
        
        function b = get.b_d(obj)
            obj.mustBeSet(obj.nw,'nw')
            obj.mustBeSet(obj.w_max,'w_max')
            b = obj.w_max*ones(2*obj.nw,1);
        end

        function F = get.F_theta(obj)
            obj.mustBeSet(obj.theta_v,'theta_v')
            np = size(obj.theta_v,1);
            F = [eye(np); -eye(np)];
        end
        
        function b = get.b_theta(obj)
            obj.mustBeSet(obj.theta_v,'theta_v')
            b = [obj.theta_v(:,2); -obj.theta_v(:,1)];
        end
        
        function v = get.disturbance_v(obj)
            obj.mustBeSet(obj.nw,'nw')
            obj.mustBeSet(obj.w_max,'w_max')
            v = obj.w_max*[-ones(obj.nw,1), ones(obj.nw,1)];
        end

        function tf = get.param_uncertainty(obj)
            tf = any(obj.theta_v);
        end

        function s = struct(obj)
            s.label = obj.label;
            s.version = obj.version;
            s.T = obj.T;
            s.N = obj.N;
            s.dt = obj.dt;
            s.x0 = obj.x0;
            s.x_ref = obj.x_ref;
            s.u_ref = obj.u_ref;
            s.F_u = obj.F_u;
            s.b_u = obj.b_u;
            s.F_x = obj.F_x;
            s.b_x = obj.b_x;
            s.h_obs = obj.h_obs;
            s.F_d = obj.F_d;
            s.b_d = obj.b_d;
            s.w_max = obj.w_max;
            s.disturbance_v = obj.disturbance_v;
            s.theta_v = obj.theta_v;
            s.theta_true = obj.theta_true;
            s.param_uncertainty = obj.param_uncertainty;
            s.Q_cost = obj.Q_cost;
            s.R_cost = obj.R_cost;
            s.P_cost = obj.P_cost;
            s.regularizer = obj.regularizer;
            s.integrator = obj.integrator;
        end
    end
    
    methods (Static,Access=private,Hidden)
        function checkPolytopes(F)
            if any(~any(F < 0) | ~any(F > 0))
                id = 'Parameters:UnboundedPolytope';
                error(id,'Polytope is unbounded.')
            end
        end
        
        function mustBePositiveDefinite(A)
            try if any(A,'all'); chol(A); end
                return
            catch
                id = 'Parameters:NotPositiveDefinite';
                error(id,'Matrix must be positive definite.')
            end
        end

        function mustBeSet(value,name)
            if isnan(value)
                id = 'Parameters:UnsetParameter';
                error(id,'The parameter %s is not set. Please set this parameter.',name)
            end
        end
    end
end