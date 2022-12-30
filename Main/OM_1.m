%excercise 1 Two body problem
clc; clear all; close all;
% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
% Initial condition
r0 = [ 26578.137; 0; 0 ]; % [km]
v0 = [ 0; 2.221; 3.173 ]; % [km/s]
y0 = [ r0; v0 ];
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 2*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );
% Plot the results
time=5*tspan;
figure(1)
hold on;
earth_sphere;
plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;
h=cross(Y(:,1:3),Y(:,4:6));
r=sqrt(Y(:,1).^2+Y(:,2).^2+Y(:,3).^2);
e=1./mu_E.*cross(Y(:,4:6),h)-Y(:,1:3)./r;
habs=zeros(length(tspan),1);
e_abs=zeros(length(tspan),1);
for i=1:length(tspan)
     habs(i)=norm(h(i,:));   %fix the norm
     e_abs(i)=norm(e(i,:));
end

figure(2);
plot(time,h,'--',time,habs,'--','linewidt',2);
xlabel('time');
ylabel('Angular momentum');
legend('h_x','h_y','h_z','|h|');

figure(3);
plot(time,e,'--',time,e_abs,'--','linewidt',2);
legend('e_x','e_y','e_z','|e|');
xlabel('time');
ylabel('Eccentricity');


figure(4);
e_dot_h=zeros(length(tspan),1);
for i=1:length(tspan)
    e_dot_h(i)=dot(e(i,:),h(i,:));
end
plot(time,e_dot_h);
xlabel('Time');
ylabel('e\cdoth');
title('Perp. check');


figure(5);
for i=1:length(tspan)
    se(i)=norm(Y(i,4:6))^2/2-mu_E/norm(Y(i,1:3));
end
plot(time,se);
title('Specific energy');
xlabel('Time');
ylabel('Specific energy');


figure(6);
ur=Y(:,1:3)./r;
uh=h./habs;
ut=cross(uh,ur);
for i=1:length(tspan)
    Vr(i)=dot(Y(i,4:6),ur(i,:));
    Vt(i)=dot(Y(i,4:6),ut(i,:));
end
plot(time,Vr,time,Vt,'linewidt',2);
xlabel('time');
ylabel('velocity components');
legend('Radial velocity','Trasversal velocity');
title('Radial and transversal velocity');

%% 
%excercise 1 Two body problem
clc; clear all; close all;
% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
% Initial condition
r0 =[ 6495; -970;-3622 ];
v0 = [4.752 ;2.130; 7.950 ];
y0 = [ r0; v0 ];
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 2*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );
% Plot the results
time=5*tspan;
figure(1)
hold on;
earth_sphere;
plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;
h=cross(Y(:,1:3),Y(:,4:6));
r=sqrt(Y(:,1).^2+Y(:,2).^2+Y(:,3).^2);
e=1./mu_E.*cross(Y(:,4:6),h)-Y(:,1:3)./r;
habs=zeros(length(tspan),1);
e_abs=zeros(length(tspan),1);
for i=1:length(tspan)
     habs(i)=norm(h(i,:));   %fix the norm
     e_abs(i)=norm(e(i,:));
end

figure(2);
plot(time,h,'--',time,habs,'--','linewidt',2);
xlabel('time');
ylabel('Angular momentum');
legend('h_x','h_y','h_z','|h|');

figure(3);
plot(time,e,'--',time,e_abs,'--','linewidt',2);
legend('e_x','e_y','e_z','|e|');
xlabel('time');
ylabel('Eccentricity');


figure(4);
e_dot_h=zeros(length(tspan),1);
for i=1:length(tspan)
    e_dot_h(i)=dot(e(i,:),h(i,:));
end
plot(time,e_dot_h);
xlabel('Time');
ylabel('e\cdoth');
title('Perp. check');


figure(5);
for i=1:length(tspan)
    se(i)=norm(Y(i,4:6))^2/2-mu_E/norm(Y(i,1:3));
end
plot(time,se);
title('Specific energy');
xlabel('Time');
ylabel('Specific energy');


figure(6);
ur=Y(:,1:3)./r;
uh=h./habs;
ut=cross(uh,ur);
for i=1:length(tspan)
    Vr(i)=dot(Y(i,4:6),ur(i,:));
    Vt(i)=dot(Y(i,4:6),ut(i,:));
end
plot(time,Vr,time,Vt,'linewidt',2);
xlabel('time');
ylabel('velocity components');
legend('Radial velocity','Trasversal velocity');
title('Radial and transversal velocity');

%% excercise 2
clc; close all; clear all;

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
Re=astroConstants(23); %Earth's equatorial radius;
J2=astroConstants(9); %Second zonal harmonic;

%Initial conditions
r0 = [26578.137; 0; 0]; %Km
v0 = [0; 2.221; 3.173] ;%km/s
y0 = [ r0; v0 ];
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 2*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ Tp, Yp ] = ode113( @(t,y) ode_2bpJ(t,y,mu_E, J2,Re,tspan), tspan, y0, options );
% Plot the results
time=5*tspan;
year=linspace(0,31536000,9000);
figure(1)  %two orbit in 5 period time
hold on;
earth_sphere;
plot3( Yp(:,1), Yp(:,2), Yp(:,3), '-','linewidt',2 )  %plot of perturbated orbit
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;

[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); %Plot Non pert. orbit
plot3( Y(:,1), Y(:,2), Y(:,3), '-' ,'linewidt',2) 

%Problem solved for 1 year for perturbed orbit 
[ Tp, Yp ] = ode113( @(t,y) ode_2bpJ(t,y,mu_E, J2,Re,year), year, y0, options );
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), year, y0, options );
figure(2);
plot3( Yp(:,1), Yp(:,2), Yp(:,3),Y(:,1), Y(:,2), Y(:,3), '-' ,'linewidt',1) 
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
colorbar
hold on;
earth_sphere;
legend('Perturbated orbit(J2)','Non perturbated','Terra');

title('Two-body problem orbit with 1 year perturbation');
axis equal;
grid on;
%Specific energy for perturbated orbit
hp=cross(Yp(:,1:3),Yp(:,4:6));
r=sqrt(Yp(:,1).^2+Yp(:,2).^2+Yp(:,3).^2);
ep=1./mu_E.*cross(Yp(:,4:6),hp)-Yp(:,1:3)./r;
habsp=zeros(length(year),1);
e_absp=zeros(length(year),1);
for i=1:length(year)
     habsp(i)=norm(hp(i,:));   %fix the norm
     e_absp(i)=norm(ep(i,:));
end

figure(3);
plot(year,hp,'--',year,habsp,'--','linewidt',2);
xlabel('time');
ylabel('Angular momentum perturbated');
legend('h_x','h_y','h_z','|h|');

figure(4);
plot(year,ep,'--',year,e_absp,'--','linewidt',2);
legend('e_x','e_y','e_z','|e|');
xlabel('time');
ylabel('Eccentricity');


figure(5);
e_dot_hp=zeros(length(year),1);
for i=1:length(year)
    e_dot_hp(i)=dot(ep(i,:),hp(i,:));
end
plot(year,e_dot_hp);
xlabel('Time');
ylabel('e\cdoth');
title('Perp. check');


figure(6);
for i=1:length(year)
    sep(i)=norm(Yp(i,4:6))^2/2-mu_E/norm(Yp(i,1:3));
end
plot(year,sep);
title('Specific energy');
xlabel('Time');
ylabel('Specific energy');


figure(7);
ur=Yp(:,1:3)./r;
uh=hp./habsp;
ut=cross(uh,ur);
for i=1:length(year)
    Vr(i)=dot(Yp(i,4:6),ur(i,:));
    Vt(i)=dot(Yp(i,4:6),ut(i,:));
end
plot(year,Vr,year,Vt,'linewidt',2);
xlabel('time');
ylabel('velocity components');
legend('Radial velocity','Trasversal velocity');
title('Radial and transversal velocity');
%% excercise 2 pt.2
clc; close all; clear all;

r0 = [6495; -970; -3622];
v0 = [4.752 ;2.130 ;7.950];

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
Re=astroConstants(23); %Earth's equatorial radius;
J2=astroConstants(9); %Second zonal harmonic;
y0 = [ r0; v0 ];
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 2*T, 9000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ Tp, Yp ] = ode113( @(t,y) ode_2bpJ(t,y,mu_E, J2,Re,tspan), tspan, y0, options );
% Plot the results
time=5*tspan;
year=linspace(0,31536000,9000);
figure(1)  %two orbit in 5 period time
hold on;
earth_sphere;
plot3( Yp(:,1), Yp(:,2), Yp(:,3), '-','linewidt',2 )  %plot of perturbated orbit
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;

[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); %Plot Non pert. orbit
plot3( Y(:,1), Y(:,2), Y(:,3), '-' ,'linewidt',2) 

%Problem solved for 1 year for perturbed orbit 
[ Tp, Yp ] = ode113( @(t,y) ode_2bpJ(t,y,mu_E, J2,Re,year), year, y0, options );
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), year, y0, options );
figure(2);
plot3( Yp(:,1), Yp(:,2), Yp(:,3),Y(:,1), Y(:,2), Y(:,3), '-' ,'linewidt',1) 
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
hold on;
earth_sphere;
legend('Perturbated orbit(J2)','Non perturbated','Terra');

title('Two-body problem orbit with 1 year perturbation');
axis equal;
grid on;
%Specific energy for perturbated orbit
hp=cross(Yp(:,1:3),Yp(:,4:6));
r=sqrt(Yp(:,1).^2+Yp(:,2).^2+Yp(:,3).^2);
ep=1./mu_E.*cross(Yp(:,4:6),hp)-Yp(:,1:3)./r;
habsp=zeros(length(year),1);
e_absp=zeros(length(year),1);
for i=1:length(year)
     habsp(i)=norm(hp(i,:));   %fix the norm
     e_absp(i)=norm(ep(i,:));
end

figure(3);
plot(year,hp,'--',year,habsp,'--','linewidt',2);
xlabel('time');
ylabel('Angular momentum perturbated');
legend('h_x','h_y','h_z','|h|');

figure(4);
plot(year,ep,'--',year,e_absp,'--','linewidt',2);
legend('e_x','e_y','e_z','|e|');
xlabel('time');
ylabel('Eccentricity');


figure(5);
e_dot_hp=zeros(length(year),1);
for i=1:length(year)
    e_dot_hp(i)=dot(ep(i,:),hp(i,:));
end
plot(year,e_dot_hp);
xlabel('Time');
ylabel('e\cdoth');
title('Perp. check');


figure(6);
for i=1:length(year)
    sep(i)=norm(Yp(i,4:6))^2/2-mu_E/norm(Yp(i,1:3));
end
plot(year,sep);
title('Specific energy');
xlabel('Time');
ylabel('Specific energy');


figure(7);
ur=Yp(:,1:3)./r;
uh=hp./habsp;
ut=cross(uh,ur);
for i=1:length(year)
    Vr(i)=dot(Yp(i,4:6),ur(i,:));
    Vt(i)=dot(Yp(i,4:6),ut(i,:));
end
plot(year,Vr,year,Vt,'linewidt',2);
xlabel('time');
ylabel('velocity components');
legend('Radial velocity','Trasversal velocity');
title('Radial and transversal velocity');
%% Excercise 3
clc; close all;clear all;
a=7000;
E0=0;
t0=0;
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
k=2;
N=100;
e=[0,0.2,0.4,0.6,0.8,0.95];
T=2*pi*sqrt(a^3/mu_E);
t=linspace(0,k*T,N);
[E]=Ecc_anom_zero(t,e,a,mu_E,t0);
figure(1)
plot(t/T,E*180/pi) 
title('Ecc. anomaly variation (Zero)');
xlabel('Period');
ylabel('E[deg]');   
legend('e=0','e=0.2,','e=0.4','e=0.6','e=0.8','e=0.95'); %plot of E(deg) respect to orbit period T
figure(2)
[E]=Ecc_anom_solve(t,e,a,mu_E,t0);
plot(t/T,E*180/pi)
title('Ecc. anomaly variation (solve)');
xlabel('Period');
ylabel('E[deg]');  
legend('e=0','e=0.2,','e=0.4','e=0.6','e=0.8','e=0.95');
figure(3)
surf(t/T,e,E'*180/pi);
