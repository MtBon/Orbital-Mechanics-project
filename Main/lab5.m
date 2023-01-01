% Lab 5
close all; clearvars; clc;
a = 15710;
e = 0.01;
i = deg2rad(87.9);
Om = pi;
om = pi;
theta = 0;
mu = astroConstants(13);
Re = astroConstants(23);
mu_moon = astroConstants(20);
J2 = astroConstants(9);
initial = [2023 1 7 10 0 0];

T = 2 * pi * sqrt(a^3/mu);
time = linspace(0,10 * T, 5000);
[r0,v0] = kp2rv(a,e , i, Om,om,theta,mu);
s0 = [r0, v0];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
%cart
[t,S] = ode113( @(t,y) ode_2bpJ(t,y,mu,J2,Re,initial,mu_moon), time, s0, options );

kep0 = [a,e,i,Om,om,theta];
%gauss
[ ~ , S_g] = ode113( @(t, s) eq_motions(t, s,J2,Re,mu,mu_moon,initial),time,kep0,options);

for u = 1 : length(time)
    [a_c(u), e_c(u), i_c(u), Om_c(u), om_c(u),theta_c(u)] = rv2kp(S(u,1:3),S(u,4:6),mu);
end
Om_c = unwrap(Om_c);
om_c = unwrap(om_c);
theta_c = unwrap(theta_c);

figure(1)
subplot(2,1,1)
plot(time/T,S_g(:,1));
subplot(2,1,2)
plot(time/T,a_c);

figure(2)
subplot(2,1,1)
plot(time/T,S_g(:,2));
subplot(2,1,2)
plot(time/T,e_c);

figure(3)
subplot(2,1,1)
plot(time/T,rad2deg(S_g(:,3)));
subplot(2,1,2)
plot(time/T,rad2deg(i_c));

figure(4)
subplot(2,1,1)
plot(time/T,rad2deg(S_g(:,4)));
subplot(2,1,2)
plot(time/T,rad2deg(Om_c));

figure(5)
subplot(2,1,1)
plot(time/T,rad2deg(S_g(:,5)));
subplot(2,1,2)
plot(time/T,rad2deg(om_c));

figure(6)
subplot(2,1,1)
plot(time/T,rad2deg(S_g(:,6)));
subplot(2,1,2)
plot(time/T,rad2deg(theta_c));


