function ds = eq_motions(t, s, J2,Re,mu,mu_moon,initial)
%
%
% Funtion to compute the derivative states of the Keplerian parameters
% using RSW  Reference System:
% INPUT : 
%           - t : time
%           - s : state [6x1]
%           - J2 : Gravitatonal field constant of the Earth 
%           - Re : Earth radius [Km]
%           - mu : Planetary constants of Earth [km^3/s^2]
%           - mu_moon : Planetary constants of Moon [km^3/s^2]
%           - initial : [1x6] Date of Departure
%
% OUTPUT :
%
%           - ds = [ d_a
%                    d_e
%                    d_i
%                    d_Omega
%                    d_omega                                                            
%                    d_theta ]
%
% CALLED FUNCTIONS : 
%                       - kp2rv ;
%                       - perturbations_fun ;
%
% AUTHORS : Matteo Bono

a = s(1);
e = s(2);
i = s(3);
Omega = s(4);
omega = s(5);
theta = s(6);

 [r , ~ ] = kp2rv(a,e,i,Omega,omega,theta,mu);

 a_pert = perturbations_fun(t,r,J2,Re,initial,mu, mu_moon);
 
%Pass in RSW reference frame for Gauss Equations
A = [ -sin(Omega) * cos(i) * sin(theta + omega) + cos(Omega) * cos(theta + omega) , cos(Omega) * cos(i) * sin(omega + theta) + sin(Omega) * cos(theta + omega) , sin(i) * sin(theta + omega)
      -sin(Omega) * cos(i) * cos(theta + omega) - cos(Omega) * sin(theta + omega) , cos(Omega) * cos(i) * cos(omega + theta) - sin(Omega) * sin(theta + omega) , sin(i) * cos(theta + omega)
      sin(Omega) * sin(i) , -cos(Omega) * sin(i) , cos(i)] ;


a_rsw = A * a_pert;

a_r = a_rsw(1);
a_s = a_rsw(2);
a_w = a_rsw(3);

 p = a * (1 - e^2);
 h = sqrt(p * mu);
 hnorm = h;
 rnorm = norm(r);
  

%Set the derivative of the Keplerian elements

     ds = [   2 * a^2/h * (e * sin(theta) * a_r + p/rnorm * a_s) 

              1/h * (p * sin(theta) * a_r + ((p + rnorm)* cos(theta) + rnorm * e) * a_s) 
            
              (rnorm * cos(theta + omega))/ h * a_w
            
              (rnorm * sin(theta + omega))/ (h * sin(i)) * a_w 
            
              1/(h * e) * (-p * cos(theta) * a_r + (p + rnorm) * sin(theta) * a_s) - (rnorm * sin(theta + omega) * cos(i) * a_w)/(hnorm * sin(i)) 
            
               h/rnorm^2 + 1/(e * hnorm) * (p * cos(theta) * a_r - (p + rnorm) * sin(theta) * a_s) ] ;
             
   

                                                                                 