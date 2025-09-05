function xdot = S1_FUAV(x, u, env, params)
% S1_FUAV - Corrected nonlinear UAV flight dynamics model (12/13-state)
%
% INPUTS
%   x - state vector [VT; alpha; beta; phi; theta; psi; P; Q; R; Pn; Pe; h]
%       angles in degrees, rates in deg/s, positions in m, VT in m/s
%   u - control vector [de; da; dr; throttle] (degrees, degrees, degrees, 0-1)
%   env - environment struct (optional)
%       .alt - ambient altitude (m)
%       .wind_body - wind vector in body frame [u; v; w] (m/s)
%   params - parameters struct
%       .aero  - aerodynamic tables (required)
%       .model - physical parameters (optional)
%       .verbose - logical
%
% OUTPUT
%   xdot - state derivative vector (12 or 13 elements as input)

% --- Input checks
if nargin < 4 || isempty(params)
    error('S1_FUAV:MissingParams', 'params structure with .aero is required');
end
if ~isfield(params,'aero') || isempty(params.aero)
    error('S1_FUAV:MissingAero', 'params.aero structure is required');
end

% Constants
DEG2RAD = pi/180;
RAD2DEG = 180/pi;
EPS = 1e-9;

% Unpack state
nstate = numel(x);
if ~ismember(nstate, [12,13])
    error('S1_FUAV:InvalidState', 'State vector x must be length 12 or 13. Received %d.', nstate);
end
VT = x(1);        % m/s
alpha_deg = x(2); beta_deg = x(3);
phi_deg = x(4); theta_deg = x(5); psi_deg = x(6);
P_deg = x(7); Q_deg = x(8); R_deg = x(9);
Pn = x(10); Pe = x(11); h = x(12);
if nstate == 13
    Thrust_state = x(13);
else
    Thrust_state = [];
end

% Controls
de = u(1); da = u(2); dr = u(3); thr = u(4);

% Environment
if isfield(env,'alt') && ~isempty(env.alt)
    alt_amb = env.alt;
else
    alt_amb = h;
end
if isfield(env,'wind_body') && numel(env.wind_body)==3
    wind_body = env.wind_body(:);
else
    wind_body = [0;0;0];
end

% Model parameters (defaults)
if isfield(params,'model') && ~isempty(params.model)
    M = params.model;
else
    M.mass = 540;
    M.S = 11.42;
    M.AR = 9.8;
    M.b = sqrt(M.AR * M.S);
    M.c = M.S / M.b;
    M.Ixx = 753.0435833333332;
    M.Iyy = 757.8323488976016;
    M.Izz = 1491.1656822309349;
    M.Ixz = -2.21398;
    M.T0 = 288.15; M.P0 = 101325; M.L = 0.0065; M.R = 287.05287; M.gamma = 1.4; M.g0 = 9.80665;
    M.engine_time_constant = 0.2;
end

% Convert to radians where needed
P = P_deg * DEG2RAD; Q = Q_deg * DEG2RAD; R = R_deg * DEG2RAD;
phi = phi_deg * DEG2RAD; theta = theta_deg * DEG2RAD; psi = psi_deg * DEG2RAD;
alpha = alpha_deg * DEG2RAD; beta = beta_deg * DEG2RAD;

% Atmosphere
[rho, ~, a, ~] = local_atmosphere(alt_amb, M);
qbar = 0.5 * rho * VT^2;
Mach = VT / a;

% qhat (nondim pitch rate)
qhat = Q * M.c / (2 * max(VT,EPS));

% Aerodynamic coefficients - use robust getter
try
    CL = get_coeff('CL', alpha_deg, de, params.aero);
catch ME
    warning('S1_FUAV:GetCL','%s', ME.message);
    CL = 0.0;
end
try
    % CD may be function of alpha & Mach or CL & Mach
    if isfield(params.aero, 'CD') && (isfield(params.aero.CD,'Mach_range') || isfield(params.aero.CD,'Mach'))
        CD = get_coeff('CD', alpha_deg, Mach, params.aero);
    else
        CL_curr = CL;
        CD = get_coeff('CD', CL_curr, Mach, params.aero);
    end
catch ME
    warning('S1_FUAV:GetCD','%s', ME.message);
    CD = 0.02;
end
try
    CY = get_coeff('CY', beta_deg, dr, params.aero);
catch ME
    warning('S1_FUAV:GetCY','%s', ME.message);
    CY = 0.0;
end
try
    Cl = get_coeff('Cl', beta_deg, da, params.aero);
catch ME
    warning('S1_FUAV:GetCl','%s', ME.message);
    Cl = 0.0;
end
try
    Cm = get_coeff('Cm', alpha_deg, qhat, de, params.aero);
catch ME
    warning('S1_FUAV:GetCm','%s', ME.message);
    Cm = 0.0;
end
try
    Cn = get_coeff('Cn', beta_deg, dr, da, params.aero);
catch ME
    warning('S1_FUAV:GetCn','%s', ME.message);
    Cn = 0.0;
end
try
    Thrust = get_coeff('Thrust', alt_amb, Mach, thr, params.aero);
catch ME
    warning('S1_FUAV:GetThrust','%s', ME.message);
    Thrust = 0;
end
Thrust = max(0, Thrust);

% Body velocity components (including wind in body frame)
u_b = VT * cos(alpha) * cos(beta) + wind_body(1);
v_b = VT * sin(beta) + wind_body(2);
w_b = VT * sin(alpha) * cos(beta) + wind_body(3);

% Recompute true airspeed from body components
VT_true = sqrt(u_b^2 + v_b^2 + w_b^2);
if VT_true < EPS, VT_true = EPS; end

% Recompute aerodynamic angles from body velocities (override inputs for consistency)
alpha = atan2(w_b, u_b);
beta = asin( max(-1,min(1, v_b / VT_true)) );
alpha_deg = alpha * RAD2DEG; beta_deg = beta * RAD2DEG;

% Aerodynamic forces (wind-axis)
Lift = CL * qbar * M.S;
Drag = CD * qbar * M.S;
Side = CY * qbar * M.S;

% Rotation from wind to body (use radians)
ca = cos(alpha); sa = sin(alpha);
cb = cos(beta); sb = sin(beta);
Rz_beta = [cb, -sb, 0; sb, cb, 0; 0, 0, 1];
Ry_alpha = [ca, 0, sa; 0, 1, 0; -sa, 0, ca];
R_body_to_wind = Rz_beta * Ry_alpha; % transforms body->wind

% Wind-frame forces (wind axes: X_w along velocity, Z_w positive downward?)
F_wind = R_body_to_wind * [Lift; Side; -Drag];
X_wind = F_wind(1); Y_wind = F_wind(2); Z_wind = F_wind(3);

% Body-frame forces (transform back)
F_body = R_body_to_wind' * F_wind;
X_body = F_body(1); Y_body = F_body(2); Z_body = F_body(3);

% Moments (body-axis)
L = Cl * qbar * M.S * M.b;
M_mom = Cm * qbar * M.S * M.c;
N = Cn * qbar * M.S * M.b;

% Inertia matrix
Ixx = M.Ixx; Iyy = M.Iyy; Izz = M.Izz; Ixz = M.Ixz;
I = [Ixx, 0, -Ixz; 0, Iyy, 0; -Ixz, 0, Izz];

% Angular velocity vector (rad/s)
omega = [P; Q; R];
cor_term = cross(omega, I*omega);
omega_dot = I \ ([L; M_mom; N] - cor_term);

% Convert angular accelerations to deg/s^2 for state storage
P_dot = omega_dot(1) * RAD2DEG;
Q_dot = omega_dot(2) * RAD2DEG;
R_dot = omega_dot(3) * RAD2DEG;

% Body accelerations approx (neglecting added mass): use Newton's 2nd law
% Add thrust along body x-axis
X_total = X_body + Thrust * cos(alpha) * cos(beta);
Y_total = Y_body + Thrust * sin(beta); % small lateral thrust comp
Z_total = Z_body - Thrust * sin(alpha) * cos(beta);

u_dot = X_total / M.mass + ( - M.g0 * sin(theta) );
v_dot = Y_total / M.mass + ( M.g0 * cos(theta) * sin(phi) );
w_dot = Z_total / M.mass + ( M.g0 * cos(theta) * cos(phi) );

% Approximate VT_dot and angle rates from body accelerations
VT_dot = (u_b * u_dot + v_b * v_dot + w_b * w_dot) / max(VT_true,EPS);

alpha_dot_rad = (w_dot * u_b - u_dot * w_b) / (max(VT_true,EPS)^2 * cos(beta));
beta_dot_rad = (v_dot * VT_true - v_b * (u_b * u_dot + v_b * v_dot + w_b * w_dot) / VT_true) / (VT_true^2);

% Add kinematic coupling terms for alpha/beta (pitch/roll rates)
alpha_dot_rad = alpha_dot_rad + Q * cos(beta) - R * sin(beta);
beta_dot_rad = beta_dot_rad + R * cos(alpha) - P * sin(alpha);

% Euler angle kinematics: compute derivatives in rad/s then convert to deg/s
phi_dot_rad = P + tan(theta) * (Q * sin(phi) + R * cos(phi));
theta_dot_rad = Q * cos(phi) - R * sin(phi);
psi_dot_rad = (Q * sin(phi) + R * cos(phi)) / max(cos(theta), 1e-6);

% Position derivatives (inertial) using rotation from body to inertial
% Rotation matrix body->inertial from phi,theta,psi (rad)
cphi = cos(phi); sphi = sin(phi);
ctheta = cos(theta); stheta = sin(theta);
cpsi = cos(psi); spsi = sin(psi);
R_bi = [ctheta*cpsi, sphi*stheta*cpsi - cphi*spsi, cphi*stheta*cpsi + sphi*spsi; ...
        ctheta*spsi, sphi*stheta*spsi + cphi*cpsi, cphi*stheta*spsi - sphi*cpsi; ...
        -stheta,     sphi*ctheta,                 cphi*ctheta];
V_body = [u_b; v_b; w_b];
V_inertial = R_bi * V_body;
Pn_dot = V_inertial(1); Pe_dot = V_inertial(2); h_dot = -V_inertial(3);

% Assemble xdot with consistent units
xdot = zeros(nstate,1);
xdot(1) = VT_dot;
xdot(2) = alpha_dot_rad * RAD2DEG;
xdot(3) = beta_dot_rad * RAD2DEG;
xdot(4) = phi_dot_rad * RAD2DEG;
xdot(5) = theta_dot_rad * RAD2DEG;
xdot(6) = psi_dot_rad * RAD2DEG;
xdot(7) = P_dot;
xdot(8) = Q_dot;
xdot(9) = R_dot;
xdot(10) = Pn_dot;
xdot(11) = Pe_dot;
xdot(12) = h_dot;

% Thrust state dynamics (if present)
if nstate == 13
    xdot(13) = (Thrust - Thrust_state) / M.engine_time_constant;
end

% Debug prints
if isfield(params,'verbose') && params.verbose
    fprintf('S1_FUAV: VT=%.2f, alpha=%.2f, beta=%.2f, CL=%.4f, CD=%.4f, Cm=%.4f\n', VT_true, alpha_deg, beta_deg, CL, CD, Cm);
end
end

% ---------------- Helper: get_coeff ----------------
function coeff = get_coeff(name, varargin)
% see previous project implementation - robust retrieval from aero struct
% Last arg must be aero struct
aero = varargin{end};
query = varargin(1:end-1);
if ~isfield(aero, name)
    error('GET_COEFF:MissingField','Aero struct does not contain %s', name);
end
A = aero.(name);
% prefer eval
if isfield(A,'eval') && isa(A.eval,'function_handle')
    try
        coeff = A.eval(query{:}); return;
    catch
        % fallback
    end
end
% interpolant
if isfield(A,'interp') && isa(A.interp,'griddedInterpolant')
    try
        coeff = A.interp(query{:}); return;
    catch
        % fallback
    end
end
% fallback to table + interpn
if isfield(A,'table') || isfield(A,'table_steady')
    if isfield(A,'table'), T = A.table; else T = A.table_steady; end
    nin = numel(query);
    switch nin
        case 2
            % alpha,de or alpha,Mach or beta,dr etc
            if isfield(A,'alpha_range') && isfield(A,'de_range')
                g1 = A.alpha_range; g2 = A.de_range;
            elseif isfield(A,'alpha_range') && isfield(A,'Mach_range')
                g1 = A.alpha_range; g2 = A.Mach_range;
            elseif isfield(A,'beta_range') && isfield(A,'dr_range')
                g1 = A.beta_range; g2 = A.dr_range;
            else
                error('GET_COEFF:GridsUnknown','Cannot determine grids for %s', name);
            end
            coeff = interpn(g1,g2,T,query{1},query{2}, 'linear', NaN);
        case 3
            if isfield(A,'alpha_range') && isfield(A,'qhat_range') && isfield(A,'de_range')
                g1=A.alpha_range; g2=A.qhat_range; g3=A.de_range;
            elseif isfield(A,'altitude_range') && isfield(A,'Mach_range') && isfield(A,'throttle_range')
                g1=A.altitude_range; g2=A.Mach_range; g3=A.throttle_range;
            elseif isfield(A,'alt_range') && isfield(A,'Mach_range') && isfield(A,'thr_range')
                g1=A.alt_range; g2=A.Mach_range; g3=A.thr_range;
            else
                error('GET_COEFF:GridsUnknown','Cannot determine 3D grids for %s', name);
            end
            coeff = interpn(g1,g2,g3,T,query{1},query{2},query{3}, 'linear', NaN);
        otherwise
            error('GET_COEFF:UnsupportedInputs','Unsupported number of inputs (%d) for %s', nin, name);
    end
    return;
end
error('GET_COEFF:NoMethod','Could not retrieve coefficient %s', name);
end

% ---------------- Helper: local_atmosphere ----------------
function [rho, P, a, T] = local_atmosphere(h, M)
if nargin < 2, M = struct(); end
if ~isfield(M,'T0'), M.T0 = 288.15; end
if ~isfield(M,'P0'), M.P0 = 101325; end
if ~isfield(M,'L'), M.L = 0.0065; end
if ~isfield(M,'R'), M.R = 287.05287; end
if ~isfield(M,'g0'), M.g0 = 9.80665; end
if ~isfield(M,'gamma'), M.gamma = 1.4; end
if h <= 11000
    T = M.T0 - M.L * h;
    P = M.P0 * (T / M.T0)^(M.g0/(M.R*M.L));
else
    T11 = 216.65; P11 = M.P0 * (T11 / M.T0)^(M.g0/(M.R*M.L));
    T = T11; P = P11 * exp(-M.g0*(h-11000)/(M.R*T11));
end
rho = P / (M.R * T);
a = sqrt(M.gamma * M.R * T);
end