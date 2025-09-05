function [x_trim, u_trim, success, resnorm, info] = solve_trim_full_s1(target, params, options)
% SOLVE_TRIM_FULL_S1 - Multi-axis trim solver using S1_FUAV dynamics
% 
% USAGE:
%   [x_trim,u_trim,success,resnorm,info] = solve_trim_full_s1(target, params, options)
% 
% INPUTS:
%   target: struct with desired flight targets (fields optional)
%     .VT    - target airspeed (m/s)
%     .h     - target altitude (m)
%     .gamma - flight-path angle (deg) positive = climb
%     .phi   - bank angle (deg) for coordinated turns
%   params: struct required by S1_FUAV with fields:
%     .aero  - aerodynamic tables (used by S1_FUAV via get_coeff)
%     .model - physical model (mass, S, c, b, inertia, g0, etc.)
%   options: optional struct with solver settings:
%     .solver = 'lsqnonlin' (default) or 'fsolve'
%     .tol = 1e-8
%     .maxiter = 200
%     .verbose = false
% 
% OUTPUTS:
%   x_trim - trimmed state vector (12 or 13 states as used by S1_FUAV)
%   u_trim - trimmed controls [de; da; dr; throttle]
%   success - logical true if solver reports success or small residual
%   resnorm - 2-norm of residual at solution
%   info - struct with solver output (exitflag, output, residual)
% 
% NOTES:
% - This solver forms a nonlinear least-squares problem by selecting a set
%   of dynamic equations from S1_FUAV and (optionally) target constraints.
% - The implementation uses a practical unknown vector consisting of
%   the primary aerodynamic and control variables. It's designed to be
%   robust and to provide a usable trim in most cases; however, you may
%   wish to tune bounds and initial guesses for specific aircraft.
% 
% Author: Generated and integrated into your repo by assistant
% Date: 2025-09-05

if nargin < 2 || isempty(params)
    error('solve_trim_full_s1:MissingParams', 'params struct with .aero and optionally .model is required');
end
if nargin < 3, options = struct(); end
if ~isfield(options,'solver'), options.solver = 'lsqnonlin'; end
if ~isfield(options,'tol'), options.tol = 1e-8; end
if ~isfield(options,'maxiter'), options.maxiter = 200; end
if ~isfield(options,'verbose'), options.verbose = false; end

% default target values
if ~isfield(target,'VT'), target.VT = 40; end
if ~isfield(target,'h'),  target.h  = 3000; end
if ~isfield(target,'gamma'), target.gamma = 0; end
if ~isfield(target,'phi'), target.phi = 0; end

% prepare model defaults if missing
if ~isfield(params,'model') || isempty(params.model)
    M = struct(); M.mass = 540; M.S = 11.42; M.c = 1.2; M.b = 10; M.g0 = 9.80665; params.model = M;
end

% Unknown vector z: [alpha_deg; beta_deg; phi_deg; theta_deg; P_deg; Q_deg; R_deg; de; da; dr; thr]
nvar = 11;

% initial guess
alpha0 = 5; beta0 = 0; phi0 = target.phi; theta0 = max(0, target.gamma + alpha0); % rough relation
P0 = 0; Q0 = 0; R0 = 0;
de0 = 0; da0 = 0; dr0 = 0; thr0 = 0.7;
z0 = [alpha0; beta0; phi0; theta0; P0; Q0; R0; de0; da0; dr0; thr0];

% bounds (deg and throttle 0..1)
lb = [-20; -20; -45; -30; -200; -200; -200; -40; -40; -40; 0];
ub = [ 40;  20;  45;  45;  200;  200;  200;  40;  40;  40; 1];

% residual function
function r = trim_residuals(z)
    % unpack
    alpha = z(1); beta = z(2); phi = z(3); theta = z(4);
    Pdeg = z(5); Qdeg = z(6); Rdeg = z(7);
    de = z(8); da = z(9); dr = z(10); thr = z(11);

    % Construct full state vector expected by S1_FUAV
    % x = [VT; alpha; beta; phi; theta; psi; P; Q; R; Pn; Pe; h]
    VT = target.VT;
    psi = 0; Pn = 0; Pe = 0; h = target.h;
    x = zeros(12,1);
    x(1) = VT; x(2) = alpha; x(3) = beta; x(4) = phi; x(5) = theta; x(6) = psi;
    x(7) = Pdeg; x(8) = Qdeg; x(9) = Rdeg; x(10) = Pn; x(11) = Pe; x(12) = h;
    u = [de; da; dr; thr];

    % call S1_FUAV
    try
        xdot = S1_FUAV(x, u, struct('alt', h), params);
    catch ME
        % if S1_FUAV errors, return large residuals so solver avoids this region
        r = ones(9,1)*1e3;
        return;
    end

    % Select residuals: use dynamic states (1:9) except we treat positions (10:12) not used
    % We'll use VT_dot, alpha_dot, beta_dot, angular accelerations P_dot,Q_dot,R_dot
    % and two trim targets: flight-path gamma and roll phi (if requested).
    r_dyn = xdot(1:3);      % VT_dot, alpha_dot, beta_dot
    r_ang = xdot(7:9);      % P_dot, Q_dot, R_dot

    % flight-path gamma constraint: gamma ~= theta - alpha (deg)
    gamma_err = (theta - alpha) - target.gamma;
    % roll angle constraint: try to enforce requested phi (for turns)
    phi_err = phi - target.phi;

    r = [r_dyn; r_ang; gamma_err; phi_err];
    % length = 3 + 3 + 1 + 1 = 8
end

% Choose solver
use_lsq = strcmpi(options.solver, 'lsqnonlin') && exist('lsqnonlin','file');
use_fsolve = exist('fsolve','file');

info = struct();
success = false; resnorm = Inf;

if use_lsq
    opts = optimoptions('lsqnonlin','Display', ternary(options.verbose,'iter','off'), 'TolFun', options.tol, 'MaxIterations', options.maxiter);
    try
        [zopt,resnorm,residual,exitflag,output] = lsqnonlin(@trim_residuals, z0, lb, ub, opts);
        info.exitflag = exitflag; info.output = output; info.residual = residual;
        success = exitflag > 0;
    catch ME
        warning('solve_trim_full_s1:lsqnonlin_failed','lsqnonlin failed: %s', ME.message);
        zopt = z0; resnorm = Inf; info.exitflag = -999; info.output = []; info.residual = [];
    end
elseif use_fsolve
    opts = optimoptions('fsolve','Display', ternary(options.verbose,'iter','off'), 'TolFun', options.tol, 'MaxIter', options.maxiter);
    try
        [zopt, fval, exitflag, output] = fsolve(@trim_residuals, z0, opts);
        info.exitflag = exitflag; info.output = output; info.residual = fval;
        resnorm = norm(fval);
        success = exitflag > 0;
    catch ME
        warning('solve_trim_full_s1:fsolve_failed','fsolve failed: %s', ME.message);
        zopt = z0; resnorm = Inf; info.exitflag = -999; info.output = []; info.residual = [];
    end
else
    error('solve_trim_full_s1:NoSolver', 'Neither lsqnonlin nor fsolve is available. Install Optimization Toolbox or provide a solver.');
end

% Unpack solution
alpha = zopt(1); beta = zopt(2); phi = zopt(3); theta = zopt(4);
Pdeg = zopt(5); Qdeg = zopt(6); Rdeg = zopt(7);
de = zopt(8); da = zopt(9); dr = zopt(10); thr = zopt(11);

x_trim = zeros(12,1);
x_trim(1) = target.VT; x_trim(2) = alpha; x_trim(3) = beta; x_trim(4) = phi; x_trim(5) = theta;
x_trim(6) = 0; x_trim(7) = Pdeg; x_trim(8) = Qdeg; x_trim(9) = Rdeg; x_trim(10) = 0; x_trim(11) = 0; x_trim(12) = target.h;

u_trim = [de; da; dr; thr];
resnorm = norm(trim_residuals(zopt));

% final success heuristics: residual norm small
success = success || (resnorm < 1e-6);

info.final_z = zopt; info.final_resnorm = resnorm;

% Return controls as u_trim
u_trim = u_trim = u_trim = u_trim; %#ok<NASGU>
% correct variable naming for return
u_trim = [de; da; dr; thr];

end

% Ternary helper
function out = ternary(cond,a,b)
    if cond, out = a; else out = b; end
end
