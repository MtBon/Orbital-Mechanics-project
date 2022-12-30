function ds = eq_motion(t, s, J2,Re,mu,mu_moon,initial)
% 
% ds = [ d_a
%        d_e
%        d_i
%        d_Omega
%        d_omega                                                            
%        d_theta ]

a = s(1);
e = s(2);
i = s(3);
Omega = s(4);
omega = s(5);
theta = s(6);

[r,v] = kp2rv(a,e,i,Omega,omega,theta,mu); 
 a_pert = perturbations_fun(t,r,J2,Re,initial,mu, mu_moon);
 
%Change coordinate in RSW

R = (r/norm(r))';
W = (cross(r,v)/norm(cross(r,v)))';
S = (cross(W,R));

a_rsw = [ R, S, W ]' * a_pert;

a_r = a_rsw(1);
a_s = a_rsw(2);
a_w = a_rsw(3);

p = a * (1 - e^2);
hnorm = sqrt(p * mu);
rnorm = norm(r);

     ds = [ 2 * a^2/hnorm * (e * sin(theta) * a_r + p/rnorm * a_s) 

      1/hnorm * (p * sin(theta) * a_r + ((p + rnorm)* cos(theta) + rnorm * e) * a_s) 
    
      (rnorm * cos(theta + omega))/ hnorm * a_w
    
      (rnorm * sin(theta + omega))/ (hnorm * sin(i)) * a_w 
    
      1/(hnorm * e) * (-p * cos(theta) * a_r + (p + rnorm) * sin(theta) * a_s) - (rnorm * sin(theta + omega) * cos(i) * a_w)/(hnorm * sin(i)) 
    
       hnorm/rnorm^2 + 1/(e * hnorm) * (p * cos(theta) * a_r - (p + rnorm) * sin(theta) * a_s)] ;
     
  

                                                                                 