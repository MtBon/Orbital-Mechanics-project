function [alpha,delta,Lon,Lat] = GroundTrack(r,Thetag0,tspan,we,t0)

%Define Angles
for i = 1:length(tspan)
delta(i,1) = asin(r(i,3)/norm(r(i,:)));
alpha(i,1) = atan2(r(i,2),r(i,1));

%Theta Greenwich
Theta.g(i,1) = Thetag0 + we*(tspan(i)-t0); 

%Longitude
Lon(i,1) = alpha(i,1) - Theta.g(i,1);

%Latitude
Lat(i,1) = delta(i,1);
end

Lat = rad2deg(Lat);
Lon = wrapToPi(Lon);
Lon = rad2deg(Lon);