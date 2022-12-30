function [r,v]=kp2rv(a,e,i,omega,w,theta,mu)
p=a*(1-e^2);
r=p/(1+e*cos(theta));
rpf=r*[cos(theta);sin(theta);0];
v=sqrt(mu/p)*[-sin(theta); e+cos(theta);0];
Rw=[cos(w) sin(w) 0;-sin(w) cos(w) 0; 0 0 1];
Ri=[1 0 0;0 cos(i) sin(i); 0 -sin(i) cos(i)];
Romega=[cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0;0 0 1];
Rge2pf=Rw*Ri*Romega;
Rpf2ge=inv(Rge2pf);
r=Rpf2ge*rpf;
v=Rpf2ge*v;
r = r';
v = v';
