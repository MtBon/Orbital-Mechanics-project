function a_pert = perturbations_fun(t,r,J2,Re,initial,mu, mu_moon)
%Compute the perturbations given the state of position
%  INPUT : 
%             - t : time [1]
%             - r : position vector [3x1] ;
%             - J2 : Gravitatonal field constant of the Earth ;
%             - Re : Earth radius [Km] ;
%             - mu : Planetary constants of Earth [km^3/s^2] ;
%             - mu_moon : Planetary constants of Moon [km^3/s^2] ;
%             - initial : [1x6] Date of Departure ;
% 
%  OUTPUT :
%             -a_pert : Acceleration of perturbations (J2 + Moon) in Cartesian Coordinates [3x1]



%Compute the position of the Moon in this time step
days = date2mjd2000(initial) + s2days(t);
[x_moon, ~ ] = ephMoon(days);

rnorm = norm(r);

%J2 Perturbation 
aj2=3/2*(J2*mu*Re^2)/rnorm^4 .* [r(1)/rnorm*(5*r(3).^2/rnorm^2-1); ...
                                r(2)/rnorm*(5*r(3)^2/rnorm^2-1) ;...
                                r(3)/rnorm*(5*r(3)^2/norm(r)^2-3)]; 


%Positions vector for the Third Body Perturbation
r_sc_m = x_moon' - r;
r_e_m = x_moon';

%Moon Perturbation
a_moon = mu_moon .* (r_sc_m / norm(r_sc_m)^3 - r_e_m / norm(r_e_m)^3);


%Overall Perturbation
a_pert = aj2 + a_moon;
end
