clc; close all; clearvars;
%satellite number :29717
i=100.3151*pi/180;
Om=101.1456*pi/180;
e=0.0075717;
w=93.9313*pi/180;
theta0=269.7867*pi/180;
N=16;
mu=astroConstants(13);
n=16.14788440*(1/(3600*24)); % Mean motion [rev/s] 
% Revolution number at epoch 10674
T=2*pi/n;
a=(mu*(T/(2*pi))^2)^(1/3);
we=15.04*(pi/180)*(1/3600);
t0=0;
m=1; k=10;

nr=(k*we)/m;
Tr=2*pi/nr;
ar=(mu*(Tr/(2*pi))^2)^(1/3);
[r,v]=kp2rv(a,e,i,Om,w,theta0,mu);
[rr,vr]=kp2rv(ar,e,i,Om,w,theta0,mu);
Orb=N*T;
Orbr=N*Tr;

y0=[r;v];
y0r=[rr;vr];

tspan = linspace( 0, Orb, 20000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu), tspan, y0, options );
[ Tr, Yr ] = ode113( @(t,y) ode_2bp(t,y,mu), tspan, y0r, options );
r=Y(:,1:3);
v=Y(:,4:6);
rr=Yr(:,1:3);
vr=Yr(:,4:6);

[alpha,delta,Lon,Lat]=GroundTrack(r,theta0,tspan,we,mu,t0);
[alphar,deltar,Lonr,Latr]=GroundTrack(rr,theta0,tspan,we,mu,t0);

Latr=Latr*180/pi;
Lonr=wrapToPi(Lonr)*180/pi;

Lat=Lat*180/pi;
Lon=wrapToPi(Lon)*180/pi;
    
figure(1)
EarthGT = imread('Earth.png');
imagesc([-180 180], [-90 90] ,flip(flip(EarthGT)));
grid on; 
axis on;
hold on;
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
set(gca,'xtick',(-180:30:180));
set(gca,'ytick',(-90:30:90))
  plot(Lon,Lat,'.b'); 
 plot(Lonr,Latr,'.r');
 plot(Lon(1),Lat(1),'dy',Lon(end),Lat(end),'^y','linewidth',8);
   plot(Lonr(1),Latr(1),'dg',Lonr(end),Latr(end),'^g','linewidth',4);
  
legend('Ground Track','Repeating','Start','End','Start repeat','End repeat');

figure(2)
earth_sphere;
hold on;
 DisegnaOrbita(a,e,w,i,Om,mu);
hold on;
  DisegnaOrbita(ar,e,w,i,Om,mu);
