function dy = ode_2bpJ( t, y, mu,J2, Re, initial,mu_moon)
%
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
% J2  Second zonal harmonic
% Re Equatorial radius
% initial Date of departure
% mu_moon Moon planetary constant

% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% CONTRIBUTORS:
% Juan Luis Gonzalo Gomez
%
% VERSIONS
% 2018-09-26: First version
% 2022-12-28: Second version including Moon perturbation
% -------------------------------------------------------------------------

% Position and velocity
r = y(1:3);
v = y(4:6);


rnorm = norm(r);

% Add the perturbation term from the oblateness of the Earth and from Moon
a_per = perturbations_fun(t,r,J2,Re,initial,mu,mu_moon);

% Set the derivatives of the state
dy = [ v
(-mu/rnorm^3)*r + a_per];
end