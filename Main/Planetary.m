% MAIN FOR THE PLANETARY MISSION
clearvars; close all; clc;

%Called functions:

%  -kp2rv ;
%  -rv2kp ;
%  -ephMoon ;
%  -GroundTrack ;
%  -astroConstants ;
%  -eq_motion ;
%  -a_pert ;
%  -matlab_graphics ;
%  -DisegnaOrbita ;




%Add path fot functions
addpath('..\Functions\');

%Function for Plots
matlab_graphics;

%Datas
Parameters.a = 115617;             %Semi major axi              [Km]
Parameters.e = 0.5671;             %Eccentricity                [-]
Parameters.incl = deg2rad(9.9112); %Orbit inclination           [Rad]
Parameters.k = 2;                  % Earth rotations            [-]
Parameters.m = 9;                  % Spacecraft revolutions     [-]
Parameters.Omega = deg2rad(0);     %RAAN                        [Rad]
Parameters.omega = deg2rad(0);     %Anomaly of periapsis        [Rad]
Parameters.theta = deg2rad(0);     %True anomaly                [Rad]
Parameters.thetaG0 = deg2rad(0);   %Anomaly of Greenwich        [Rad]
Parameters.mu = astroConstants(13);
Parameters.Re = astroConstants(23);
Parameters.we = deg2rad(15.04 / 3600);  %Earth angular velocity [Rad/s]
mu_moon = astroConstants(20);
J2 = astroConstants(9);


%Times for different analysis
times.T = 2*pi*sqrt(Parameters.a^3/Parameters.mu);
times.day = 24 * 3600;
times.days_10 = 10 * times.day;
times.t0 = 0;

%Times initialization for GTs
T_span = linspace(0, times.T, 9000);
day_span = linspace(0, times.day, 9000);
days10_span = linspace(0, times.days_10, 9000);

%Date of Departure
initial = [2023 1 7 10 0 0];


%Cartesian 
[r0, v0] = kp2rv(Parameters.a, Parameters.e, Parameters.incl, Parameters.Omega, ...
                                        Parameters.omega, Parameters.theta, Parameters.mu);


[r_moon, v_moon] = ephMoon(date2mjd2000(initial));
y0_moon = [r_moon , v_moon];

y0 = [r0 v0];  %Initial State
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[~ , Y_nom_orbit ] = ode113( @(t,y) ode_2bp(t,y,Parameters.mu), T_span, y0, options );
[~ , Y_nom_day ] = ode113( @(t,y) ode_2bp(t,y,Parameters.mu), day_span, y0, options );
[~ , Y_nom_days10 ] = ode113( @(t,y) ode_2bp(t,y,Parameters.mu), days10_span, y0, options );

%Propagation of Moon Orbit
[~ , Y_moon_orbit ] = ode113( @(t,y) ode_2bp(t,y,Parameters.mu), T_span, y0_moon, options );
[~ , Y_moon_day ] = ode113( @(t,y) ode_2bp(t,y,Parameters.mu), day_span, y0_moon, options );
[~ , Y_moon_days10 ] = ode113( @(t,y) ode_2bp(t,y,Parameters.mu), days10_span, y0_moon, options );



%States
r_orb = Y_nom_orbit(:,1:3);
v_orb = Y_nom_orbit(:,4:6);

r_day = Y_nom_day(:,1:3);
v_day = Y_nom_day(:,4:6);

r_days10 = Y_nom_days10(:,1:3);
v_days10 = Y_nom_days10(:,4:6);


%Nominal GroundTrack for 1 Orbit
[~,~,Lon.nominal.orbit,Lat.nominal.orbit] = GroundTrack(r_orb,Parameters.thetaG0,...
                                             T_span,Parameters.we,times.t0);

%Nominal GroundTrack for 1 Day
[~,~,Lon.nominal.day,Lat.nominal.day] = GroundTrack(r_day,Parameters.thetaG0,...
                                             day_span,Parameters.we,times.t0);

%Nominal GroundTrack for 10 Days
[~,~,Lon.nominal.days10,Lat.nominal.days10] = GroundTrack(r_days10,Parameters.thetaG0,...
                                             days10_span,Parameters.we,times.t0);



% REPEATING GT


Parameters.n = (Parameters.k * Parameters.we) / Parameters.m;
times.T_rep = 2 * pi / Parameters.n;
Trep_span = linspace(0, times.T_rep, 6000);

Parameters.a_rep = (Parameters.mu * (times.T_rep/(2*pi))^2)^(1/3);

fprintf('Semi major axis for Repeating GT = %.2fKm\n',Parameters.a_rep);

[r0_rep, v0_rep] = kp2rv(Parameters.a_rep,Parameters.e,Parameters.incl,...
                                Parameters.Omega,Parameters.omega,Parameters.theta,Parameters.mu);

%Initial state for Repeating GT
y0_rep = [r0_rep, v0_rep];

[~ , Y_rep_orbit ] = ode113( @(t,y) ode_2bp(t,y,Parameters.mu), Trep_span, y0_rep, options );
[~ , Y_rep_day ] = ode113( @(t,y) ode_2bp(t,y,Parameters.mu), day_span, y0_rep, options );
[~ , Y_rep_days10 ] = ode113( @(t,y) ode_2bp(t,y,Parameters.mu), days10_span, y0_rep, options );

%States
r_rep_orb = Y_rep_orbit(:,1:3);
v_rep_orb = Y_nom_orbit(:,4:6);

r_rep_day = Y_rep_day(:,1:3);
v_rep_day = Y_rep_day(:,4:6);

r_rep_days10 = Y_rep_days10(:,1:3);
v_rep_days10 = Y_rep_days10(:,4:6);


%Repeated GroundTrack for 1 Orbit
[~,~,Lon.rep.orbit,Lat.rep.orbit] = GroundTrack(r_rep_orb,Parameters.thetaG0,...
                                             Trep_span,Parameters.we,times.t0);

%Repeated GroundTrack for 1 Day
[~,~,Lon.rep.day,Lat.rep.day] = GroundTrack(r_rep_day,Parameters.thetaG0,...
                                             day_span,Parameters.we,times.t0);

%Repeated GroundTrack for 10 Days
[~,~,Lon.rep.days10,Lat.rep.days10] = GroundTrack(r_rep_days10,Parameters.thetaG0,...
                                             days10_span,Parameters.we,times.t0);



[~ , Y_pert_orbit ] = ode113( @(t,y) ode_2bpJ(t,y,Parameters.mu,J2,Parameters.Re ,initial,mu_moon), T_span, y0, options );
[~ , Y_pert_day ] = ode113( @(t,y) ode_2bpJ(t,y,Parameters.mu,J2, Parameters.Re ,initial,mu_moon), day_span, y0, options );
[~ , Y_pert_days10 ] = ode113( @(t,y) ode_2bpJ(t,y,Parameters.mu,J2,Parameters.Re , initial,mu_moon), days10_span, y0, options );


[~ , Y_rep_pert_orbit ] = ode113( @(t,y) ode_2bpJ(t,y,Parameters.mu,J2,Parameters.Re , initial,mu_moon), Trep_span, y0_rep, options );
[~ , Y_rep_pert_day ] = ode113( @(t,y) ode_2bpJ(t,y,Parameters.mu,J2, Parameters.Re ,initial,mu_moon), day_span, y0_rep, options );
[~ , Y_rep_pert_days10 ] = ode113( @(t,y) ode_2bpJ(t,y,Parameters.mu,J2,Parameters.Re , initial,mu_moon), days10_span, y0_rep, options );

%States for Repeated Perturbed GT
r_rep_pert_orb = Y_rep_pert_orbit(:,1:3);
v_rep_pert_orb = Y_rep_pert_orbit(:,4:6);

r_rep_pert_day = Y_rep_pert_day(:,1:3);
v_rep_pert_day = Y_rep_pert_day(:,4:6);

r_rep_pert_days10 = Y_rep_pert_days10(:,1:3);
v_rep_pert_days10 = Y_rep_pert_days10(:,4:6);


%Repeated Perturbed GroundTrack for 1 Orbit
[~,~,Lon.rep.per.orbit,Lat.rep.per.orbit] = GroundTrack(r_rep_pert_orb,Parameters.thetaG0,...
                                             Trep_span,Parameters.we,times.t0);

%Repeated GroundTrack for 1 Day
[~,~,Lon.rep.per.day,Lat.rep.per.day] = GroundTrack(r_rep_pert_day,Parameters.thetaG0,...
                                             day_span,Parameters.we,times.t0);

%Repeated GroundTrack for 10 Days
[~,~,Lon.rep.per.days10,Lat.rep.per.days10] = GroundTrack(r_rep_pert_days10,Parameters.thetaG0,...
                                             days10_span,Parameters.we,times.t0);

%States fot Nominal Perturbed GT
r_pert_orb = Y_pert_orbit(:,1:3);
v_pert_orb = Y_pert_orbit(:,4:6);

r_pert_day = Y_pert_day(:,1:3);
v_pert_day = Y_pert_day(:,4:6);

r_pert_days10 = Y_pert_days10(:,1:3);
v_pert_days10 = Y_pert_days10(:,4:6);

%Repeated Perturbed GroundTrack for 1 Orbit
[~,~,Lon.nominal.per.orbit,Lat.nominal.per.orbit] = GroundTrack(r_pert_orb,Parameters.thetaG0,...
                                             Trep_span,Parameters.we,times.t0);

%Repeated GroundTrack for 1 Day
[~,~,Lon.nominal.per.day,Lat.nominal.per.day] = GroundTrack(r_pert_day,Parameters.thetaG0,...
                                             day_span,Parameters.we,times.t0);

%Repeated GroundTrack for 10 Days
[~,~,Lon.nominal.per.days10,Lat.nominal.per.days10] = GroundTrack(r_pert_days10,Parameters.thetaG0,...
                                             days10_span,Parameters.we,times.t0);



%Propagation with Perturbations
times.T_prop = linspace(0, 60 * times.days_10,10000 ); %600 days propagation

%Initial States

%For Cartesian
s0 = [Parameters.a Parameters.e Parameters.incl Parameters.Omega Parameters.omega Parameters.theta];

%For Gauss Equations
kep0 = [Parameters.a, Parameters.e,Parameters.incl,Parameters.Omega,Parameters.omega,Parameters.theta];

%Gauss equations
[ ~ , S_gauss] = ode113( @(t, s) eq_motion(t, s,J2,Parameters.Re,Parameters.mu,mu_moon,initial),times.T_prop,kep0,options);

%Propagation of orbit in Cartesian
[~ , S_per ] = ode113( @(t,y) ode_2bpJ(t,y,Parameters.mu,J2,Parameters.Re , initial,mu_moon), times.T_prop, y0, options );


%Motion of Moon
s_moon = zeros(3,length(times.T_prop));

for i = 1:length(times.T_prop)
    initial(3) = initial(3) + s2days(times.T_prop(i));
[s_moon(:,i), ~ ] = ephMoon(date2mjd2000(initial));
initial(3) = initial(3) - s2days(times.T_prop(i));
end


%Variation of Kep elements in Cartesian
r_per = S_per(:,1:3);
v_per = S_per(:,4:6);
a_car = zeros(length(times.T_prop),1);
e_car = a_car;
i_car = a_car;
Om_car = a_car;
om_car = a_car;
theta_car = a_car;
r_car = zeros(length(times.T_prop),3);
v_car = r_car;
err_x = zeros(length(times.T_prop),1);
err_y = err_x;
err_z = err_x;

for i = 1:length(times.T_prop)

    [a_car(i),e_car(i),i_car(i),Om_car(i),om_car(i),theta_car(i)] = rv2kp(r_per(i,:),v_per(i,:),Parameters.mu);
    [r_car(i,:) , v_car(i,:)] = kp2rv(S_gauss(i,1),S_gauss(i,2),S_gauss(i,3),S_gauss(i,4),S_gauss(i,5),S_gauss(i,6),Parameters.mu);
    err_x(i) = r_per(i,1) - r_car(i,1);
    err_y(i) = r_per(i,2) - r_car(i,2);
    err_z(i) = r_per(i,3) - r_car(i,3);
    
end

%Unwrapping
Om_car = unwrap(Om_car);
om_car = unwrap(om_car);
theta_car = unwrap(theta_car);


%Variations of Keplerian Parameters with Gauss equations
a_prop = S_gauss(:,1);
e_prop = S_gauss(:,2);
incl_prop = S_gauss(:,3);
Omega_prop = S_gauss(:,4);
omega_prop = S_gauss(:,5);
theta_prop = S_gauss(:,6);


%GT for 1 Orbit
figure(1)
EarthGT = imread('Earth.png');
hold on;
imagesc([-180 180], [-90 90] ,(flip(EarthGT))); 

xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
set(gca,'xtick',(-180:30:180));
set(gca,'ytick',(-90:30:90));
%Nominal 
plot(Lon.nominal.orbit,Lat.nominal.orbit,'.g');
plot(Lon.nominal.orbit(1),Lat.nominal.orbit(1),'dr'...
     ,Lon.nominal.orbit(end),Lat.nominal.orbit(end),'^r','linewidth',6);
%Repeated
plot(Lon.rep.orbit,Lat.rep.orbit,'.b');
plot(Lon.rep.orbit(1),Lat.rep.orbit(1),'dy'...
     ,Lon.rep.orbit(end),Lat.rep.orbit(end),'^y','linewidth',6);

%Nominal Perturbed
plot(Lon.nominal.per.orbit,Lat.nominal.per.orbit,'.w');
plot(Lon.nominal.per.orbit(1),Lat.nominal.per.orbit(1),'dc'...
     ,Lon.nominal.per.orbit(end),Lat.nominal.per.orbit(end),'^c','linewidth',6);

%Repeated Perturbed
plot(Lon.rep.per.orbit,Lat.rep.per.orbit,'.m');
plot(Lon.rep.per.orbit(1),Lat.rep.per.orbit(1),'dm'...
     ,Lon.rep.per.orbit(end),Lat.rep.per.orbit(end),'^m','linewidth',6);

legend('Nominal','Start Nominal','End Nominal','Repeated','Start Repeated','End Repeated','Nominal perturbed'...
                              ,'Start Nom Pert','End Nom Pert','Repe Pert','Start Rep Per','End Rep Pert');
xlim([-180 180]);
ylim([-90 90]);

%GT for 1 Day
figure(2)
EarthGT = imread('Earth.png');
hold on;
imagesc([-180 180], [-90 90] ,(flip(EarthGT))); 

xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
set(gca,'xtick',(-180:30:180));
set(gca,'ytick',(-90:30:90));
%Nominal 
plot(Lon.nominal.day,Lat.nominal.day,'.g');
plot(Lon.nominal.day(1),Lat.nominal.day(1),'dr'...
     ,Lon.nominal.day(end),Lat.nominal.day(end),'^r','linewidth',6);
%Repeated
plot(Lon.rep.day,Lat.rep.day,'.b');
plot(Lon.rep.day(1),Lat.rep.day(1),'dy'...
     ,Lon.rep.day(end),Lat.rep.day(end),'^y','linewidth',6);

%Nominal Perturbed
plot(Lon.nominal.per.day,Lat.nominal.per.day,'.w');
plot(Lon.nominal.per.day(1),Lat.nominal.per.day(1),'dc'...
     ,Lon.nominal.per.day(end),Lat.nominal.per.day(end),'^c','linewidth',6);
%Repeated perturbed
plot(Lon.rep.per.day,Lat.rep.per.day,'.m');
plot(Lon.rep.per.day(1),Lat.rep.per.day(1),'dm'...
     ,Lon.rep.per.day(end),Lat.rep.per.day(end),'^m','linewidth',6);
legend('Nominal','Start Nominal','End Nominal','Repeated','Start Repeated','End Repeated','Nominal perturbed'...
                              ,'Start Nom Pert','End Nom Pert','Repe Pert','Start Rep Per','End Rep Pert');

xlim([-180 180]);
ylim([-90 90]);



%GT for 10 Days
figure(3)
EarthGT = imread('Earth.png');
hold on;
imagesc([-180 180], [-90 90] ,(flip(EarthGT))); 

xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
set(gca,'xtick',(-180:30:180));
set(gca,'ytick',(-90:30:90));
%Nominal
plot(Lon.nominal.days10,Lat.nominal.days10,'.g');
plot(Lon.nominal.days10(1),Lat.nominal.days10(1),'dr'...
     ,Lon.nominal.days10(end),Lat.nominal.days10(end),'^r','linewidth',6);
%Repeated
plot(Lon.rep.days10,Lat.rep.days10,'.b');
plot(Lon.rep.days10(1),Lat.rep.days10(1),'dy'...
     ,Lon.rep.days10(end),Lat.rep.days10(end),'^y','linewidth',6);

%Nominal Perturbed
plot(Lon.nominal.per.days10,Lat.nominal.per.days10,'.w');
plot(Lon.nominal.per.days10(1),Lat.nominal.per.days10(1),'dc'...
     ,Lon.nominal.per.days10(end),Lat.nominal.per.days10(end),'^c','linewidth',6);
%Repeated perturbed
plot(Lon.rep.per.days10,Lat.rep.per.days10,'.m');
plot(Lon.rep.per.days10(1),Lat.rep.per.days10(1),'dm'...
     ,Lon.rep.per.days10(end),Lat.rep.per.days10(end),'^m','linewidth',6);
legend('Nominal','Start Nominal','End Nominal','Repeated','Start Repeated','End Repeated','Nominal perturbed'...
                              ,'Start Nom Pert','End Nom Pert','Repe Pert','Start Rep Per','End Rep Pert');

xlim([-180 180]);
ylim([-90 90]);


figure(4)
earth_sphere;
hold on;
DisegnaOrbita(Parameters.a, Parameters.e, Parameters.omega, Parameters.incl, Parameters.Omega)
DisegnaOrbita(Parameters.a_rep, Parameters.e, Parameters.omega, Parameters.incl, Parameters.Omega)
plot3(y0_moon(1),y0_moon(2),y0_moon(3),'.','Markersize',11);
grid on;
legend('Earth','Nominal Orbit','Repeated Orbit','Moon');


%Orbit representation
figure(5)
plot3(S_per(:,1), S_per(:,2), S_per(:,3));
hold on;
plot3(s_moon(1,:),s_moon(2,:),s_moon(3,:));
hold on;
earth_sphere;


%Variation of Semi major axis
figure(6)
title('Variations of Semi major axis');

subplot(2,1,1)
plot(times.T_prop/times.day, a_prop);
hold on;
plot(times.T_prop/times.day, movmean(a_prop,300));
legend('Gauss','Secular');
xlabel('Time [Days]');
ylabel('a [Km]');

subplot(2,1,2)
plot(times.T_prop/times.day, a_car);
legend('Cartesian');
xlabel('Time [Days]');
ylabel('a [Km]');


%Variation of eccentricity
figure(7)
title('Variations of eccentricity');

subplot(2,1,1)
plot(times.T_prop/times.day, e_prop);
hold on;
plot(times.T_prop/times.day, movmean(e_prop,950));
legend('Gauss','Secular');
xlabel('Time [Days]');
ylabel('e [-]');

subplot(2,1,2)
plot(times.T_prop/times.day, e_car);
legend('Cartesian');
xlabel('Time [Days]');
ylabel('e [-]');


%Variation of Inclination
figure(8)
title('Variations of Inclination');

subplot(2,1,1)
plot(times.T_prop/times.day, rad2deg(incl_prop));
hold on;
plot(times.T_prop/times.day, movmean(rad2deg(incl_prop),100));
plot(times.T_prop/times.day, movmean(rad2deg(incl_prop),900));
legend('Gauss','Long Period','Secular');
xlabel('Time [Days]');
ylabel('i [Deg]');

subplot(2,1,2)
plot(times.T_prop/times.day, rad2deg(i_car));
legend('Cartesian');
xlabel('Time [Days]');
ylabel('i [Deg]');


%Variation of RAAN
figure(9)
title('Variations of RAAN');

subplot(2,1,1)
plot(times.T_prop/times.day, rad2deg(Omega_prop));
hold on;
plot(times.T_prop/times.day, movmean(rad2deg(Omega_prop),100));
plot(times.T_prop/times.day, movmean(rad2deg(Omega_prop),500));
legend('Gauss','Long Term','Secular');
xlabel('Time [Days]');
ylabel('\Omega [Deg]');

subplot(2,1,2)
plot(times.T_prop/times.day, rad2deg(Om_car));
legend('Cartesian');
xlabel('Time [Days]');
ylabel('\Omega [Deg]');

%Variations of anomaly of periapsis
figure(10)
title('Variations of anomaly of periapsis');

subplot(2,1,1)
plot(times.T_prop/times.day, rad2deg(omega_prop));
hold on;
plot(times.T_prop/times.day, movmean(rad2deg(omega_prop),100));
plot(times.T_prop/times.day, movmean(rad2deg(omega_prop),900));
legend('Gauss','Long Term','Secular');
xlabel('Time [Days]');
ylabel('\omega [Deg]');

subplot(2,1,2)
plot(times.T_prop/times.day, rad2deg(om_car));
legend('Cartesian');
xlabel('Time [Days]');
ylabel('\omega [Deg]');


%mov_orbit(S_gauss,times.T_prop);



