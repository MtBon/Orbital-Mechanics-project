% MAIN FOR THE PLANETARY MISSION
clearvars; close all; clc;

%Called functions:






%Add path fot functions
addpath('..\Functions\');


%Datas
Parameters.a = 115617;             %Semi major axi      [Km]
Parameters.e = 0.5671;             %Eccentricity     [-]
Parameters.incl = deg2rad(9.9112); %Orbit inclination       [Rad]
Parameters.k = 2;                  % Earth rotations        [-]
Parameters.m = 9;                  % Spacecraft revolutions     [-]
Parameters.Omega = deg2rad(0);     %RAAN        [Rad]
Parameters.omega = deg2rad(0);     %Anomaly of periapsis    [Rad]
Parameters.theta = deg2rad(0);     %True anomaly        [Rad]
Parameters.thetaG0 = deg2rad(0);   %Anomaly of Greenwich    [Rad]
Parameters.mu = astroConstants(13);
Parameters.Re = astroConstants(23);
Parameters.we = deg2rad(15.04 / 3600);


%Times for different analysis
times.T = 2*pi*sqrt(Parameters.a^3/Parameters.mu);
times.day = 24 * 3600;
times.days_10 = 10 * times.day;
times.t0 = 0;

%Times initialization for GTs
T_span = linspace(0, times.T, 9000);
day_span = linspace(0, times.day, 9000);
days10_span = linspace(0, times.days_10, 9000);

%Cartesian 
[r0, v0] = kp2rv(Parameters.a, Parameters.e, Parameters.incl, Parameters.Omega, ...
                                        Parameters.omega, Parameters.theta, Parameters.mu);


y0 = [r0 v0];  %Initial State
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[~ , Y_nom_orbit ] = ode113( @(t,y) ode_2bp(t,y,Parameters.mu), T_span, y0, options );
[~ , Y_nom_day ] = ode113( @(t,y) ode_2bp(t,y,Parameters.mu), day_span, y0, options );
[~ , Y_nom_days10 ] = ode113( @(t,y) ode_2bp(t,y,Parameters.mu), days10_span, y0, options );

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

legend('Nominal','Start Nominal','End Nominal','Repeated','Start Repeated','End Repeated');
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
legend('Nominal','Start Nominal','End Nominal','Repeated','Start Repeated','End Repeated');
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
legend('Nominal','Start Nominal','End Nominal','Repeated','Start Repeated','End Repeated');
xlim([-180 180]);
ylim([-90 90]);


figure(4)
earth_sphere;
hold on;
DisegnaOrbita(Parameters.a, Parameters.e, Parameters.omega, Parameters.incl, Parameters.Omega)
DisegnaOrbita(Parameters.a_rep, Parameters.e, Parameters.omega, Parameters.incl, Parameters.Omega)



