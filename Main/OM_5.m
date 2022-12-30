clc; close all; clearvars;
%Initial Conditions
mu = astroConstants(13);
Re=astroConstants(23);
J2 = astroConstants(9);
a0 = 7571;
e0 = 0.01;
i0 = deg2rad(87.9);
Omega0 = deg2rad(87.9);
omega0 = deg2rad(180);
theta0 = deg2rad(0);
Tp = 2 * pi * sqrt(a0^3/mu);
tspan = linspace(0,100*Tp,10000);
s0 = [a0;e0;i0;Omega0;omega0;theta0]';


options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

[T, S] = ode113( @(t, s) eq_motion(t, s, @(t, s) perturbations_fun(t, s, J2, Re,mu),J2,Re,mu),tspan,s0,options);


