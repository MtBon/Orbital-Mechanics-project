%% GROUND TRACK
% Data

clc; clearvars;close all;

addpath('..\IAMS\');
addpath('..\Functions\');

a=8350;
e=0.1976;
i=60*pi/180;
Om=270*pi/180;
w=45*pi/180;
theta0=230*pi/180;
Theta.g0=0;
N=3.25;  %N orbits
mu = astroConstants(13);
Re=astroConstants(23);
we=15.04*(pi/180)*(1/3600);
T=2*pi*sqrt(a^3/mu);
Orb=N*T;
t0=0;

[p,v]=kp2rv(a,e,i,Om,w,theta0,mu);
y0=[p;v]';
tspan = linspace( 0, Orb, 3000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu), tspan, y0, options );
r=Y(:,1:3);
v=Y(:,4:6);
%Define Angles
for i=1:length(tspan)
delta(i,1)=asin(r(i,3)/norm(r(i,:)));

alpha(i,1)=atan2(r(i,2),r(i,1));
%Theta Greenwich
Theta.g(i,1)= Theta.g0 + we*(tspan(i)-t0); 

%Longitude
Lon(i,1)=alpha(i,1)-Theta.g(i,1);

%Latitude
Lat(i,1)=delta(i,1);
end



Lat=Lat*180/pi;
Lon=wrapToPi(Lon)*180/pi;
figure(1)
EarthGT = imread('Earth.png'); hold on;
imagesc([-180 180], [-90 90] ,flip(EarthGT)); 
hold on;
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
set(gca,'xtick',(-180:30:180));
set(gca,'ytick',(-90:30:90));
  plot(Lon,Lat,'.b');
  hold on;
 
 hold on;
 plot(Lon(1),Lat(1),'dy',Lon(end),Lat(end),'^y','linewidth',8); hold on;
  legend('Ground Track','Start','End');

figure(2)
earth_sphere;
hold on;
 DisegnaOrbita(a,e,w,i,Om);
hold on;
figure(3)
earth_sphere;
hold on;
plot3(r(:,1),r(:,2),r(:,3)); hold on;
plot3(p(1),p(2),p(3),'.','MarkerSize',14);

%% CASE 2 
clc; close all;clearvars;
a=42166.167;  e=0.0;  i=0*pi/180;  Om=0*pi/180;   w=0*pi/180;   theta0=20*pi/180;  Theta.g0=0; 
N=500;
mu = astroConstants(13);
Re=astroConstants(23);
we=15.04*(pi/180)*(1/3600);
T=2*pi*sqrt(a^3/mu);
Orb=N*T;

t0=0;

[p,v]=kp2rv(a,e,i,Om,w,theta0,mu);
y0=[p ,v];
tspan = linspace( 0, Orb, 10000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu), tspan, y0, options );
r=Y(:,1:3);
v=Y(:,4:6);
%Define Angles
for i=1:length(tspan)
delta(i,1)=asin(r(i,3)/norm(r(i,:)));

alpha(i,1)=atan2(r(i,2),r(i,1));
%Theta Greenwich
Theta.g(i,1)= Theta.g0 + we*(tspan(i)-t0); 

%Longitude
Lon(i,1)=alpha(i,1)-Theta.g(i,1);

%Latitude
Lat(i,1)=delta(i);
end


Lat=Lat*180/pi;
Lon=wrapToPi(Lon)*180/pi;
figure(1)
EarthGT = imread('Earth.png'); hold on;
imagesc([-180 180], [-90 90] ,(flip(EarthGT))); 
hold on;
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
set(gca,'xtick',(-180:30:180));
set(gca,'ytick',(-90:30:90));
  plot(Lon,Lat,'.b');
  hold on;
 hold on;
 plot(Lon(1),Lat(1),'dy',Lon(end),Lat(end),'^y','linewidth',8); hold on;

legend('Ground Track','Start','End');

figure(2)
earth_sphere;
hold on;
 DisegnaOrbita(a,e,w,i,Om);
hold on;




%% CASE 1 Repeating GT   giusto algoritmo

clc; close all;clearvars;
%Case 1
% a=8350;   e=0.1976;   i=60*pi/180;  Om=270*pi/180;  w=45*pi/180;  theta0=230*pi/180; Theta.g0=0*pi/180; N=15;  
% k=12; m=1;
%Datas case 2
% a=26600;  e=0.74;  i=63.4*pi/180;  Om=50*pi/180;   w=280*pi/180;   theta0=0*pi/180;  Theta.g0=0*pi/180; 
% N=30;        m=1;     %Earth's revolutions;        
% k=2;   %Satellite's revolutions;


%Data case 3       
a=7171.010;  e=0;  i=98*pi/180;  Om=0*pi/180;   w=40*pi/180;   theta0=0*pi/180;  Theta.g0=0*pi/180;
N=30;        m=1;     %Earth's revolutions;        
k=15;   %Satellit

%Case example
% a=16733.65;  e=0.1976;  i=60*pi/180;  Om=270*pi/180;  w=45*pi/180; theta0=230*pi/180;
% k=4;  m=1;  N=4; Theta.g0=0*pi/180;

mu = astroConstants(13);

Re=astroConstants(23);
we=15.04*(pi/180)*(1/3600);
t0=0;



% Repeating GT
n=(k*we)/m;
Tr=2*pi/n;
T=2*pi*sqrt(a^3/mu);
ar=(mu*(Tr/(2*pi))^2)^(1/3);
Orbr=N*Tr;
Orb=N*T;

[p,v]=kp2rv(a,e,i,Om,w,theta0,mu);
[pr,vr]=kp2rv(ar,e,i,Om,w,theta0,mu);
y0=[p v];
y0r=[pr vr];

tspanr = linspace( 0, Orbr, 40000 );
tspan = linspace( 0, Orb, 40000 );
% Set options for the ODE solver

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu), tspan, y0, options );
[ Tr, Yr ] = ode113( @(t,y) ode_2bp(t,y,mu), tspanr, y0r, options );
r=Y(:,1:3);
v=Y(:,4:6);
rr=Yr(:,1:3);
vr=Yr(:,4:6);

[alpha,delta,Lon,Lat]=GroundTrack(r,Theta.g0,tspan,we,t0);
[alphar,deltar,Lonr,Latr]=GroundTrack(rr,Theta.g0,tspanr,we,t0);



figure(1)
EarthGT = imread('Earth.png');hold on;
imagesc([-180 180], [-90 90] ,(flip(EarthGT))); 

xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
set(gca,'xtick',(-180:30:180));
set(gca,'ytick',(-90:30:90));
  plot(Lon,Lat,'.b');
  hold on;
 plot(Lonr,Latr,'.r');
 hold on;
 plot(Lon(1),Lat(1),'dy',Lon(end),Lat(end),'^y','linewidth',8); hold on;
   plot(Lonr(1),Latr(1),'dg',Lonr(end),Latr(end),'^g','linewidth',4);
  
legend('Ground Track','Repeating','Start','End','Start repeat','End repeat');
xlim([-180 180]);
ylim([-90 90]);

figure(2)
earth_sphere;
hold on;
 DisegnaOrbita(a,e,w,i,Om);
hold on;
  DisegnaOrbita(ar,e,w,i,Om);
