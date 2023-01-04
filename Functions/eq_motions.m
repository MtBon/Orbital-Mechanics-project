function ds = eq_motions(t, s, J2,Re,mu,mu_moon,initial)
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
 
%RSW reference frame
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
 vnorm = norm(v);
 

% %Chane in TNH RF
% T = v/vnorm;
% H = cross(r,v);
% N = cross(H,T);
% 
% A = [T' N' H'];
% 
% a_tnh = A' * a_pert;
% 
% a_t = a_tnh(1);
% a_n = a_tnh(2);
% a_h = a_tnh(3);
% 
% 
% 
% da = 2 * a^2 * vnorm/mu * a_t;
%     de =   1/vnorm * (2 * (e + cos(theta)) * a_t - rnorm/a * sin(theta) * a_n);
%     di =   rnorm * cos(theta + omega)/h * a_h;
%      dOm =  rnorm * sin(theta + omega)/(h * sin(i)) * a_h;
%      dom =  1/( e * vnorm) * (2 * sin(theta) * a_t + (2 * e + rnorm/a * cos(theta)) * a_n) - rnorm * sin(theta + omega) * cos(i) / (h * sin(i)) * a_h;
%       dtheta =  h/rnorm^2 - 1 / (e * vnorm) * (2 * sin(theta) * a_t + (2 * e + rnorm/a * cos(theta)) * a_n);
% 
%       ds = [da, de , di , dOm , dom , dtheta]';

     ds = [ 2 * a^2/h * (e * sin(theta) * a_r + p/rnorm * a_s) 

      1/h * (p * sin(theta) * a_r + ((p + rnorm)* cos(theta) + rnorm * e) * a_s) 
    
      (rnorm * cos(theta + omega))/ h * a_w
    
      (rnorm * sin(theta + omega))/ (h * sin(i)) * a_w 
    
      1/(h * e) * (-p * cos(theta) * a_r + (p + rnorm) * sin(theta) * a_s) - (rnorm * sin(theta + omega) * cos(i) * a_w)/(hnorm * sin(i)) 
    
       h/rnorm^2 + 1/(e * hnorm) * (p * cos(theta) * a_r - (p + rnorm) * sin(theta) * a_s)] ;
     
%   

                                                                                 