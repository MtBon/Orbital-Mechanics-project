function a_pert = perturbations_fun(t,r,J2,Re,initial,mu, mu_moon)

%Compute the perturbations given the state
initial(3) = initial(3) + s2days(t);
[x_moon, ~ ] = ephMoon(date2mjd2000(initial));
initial(3) = initial(3) - s2days(t);

rnorm = norm(r);

aj2=3/2*(J2*mu*Re^2)/rnorm^4 .* [r(1)/rnorm*(5*r(3).^2/rnorm^2-1); ...
                                r(2)/rnorm*(5*r(3)^2/rnorm^2-1) ;...
                                r(3)/rnorm*(5*r(3)^2/norm(r)^2-3)]; 

r_sc_m = x_moon' - r;
r_e_m = x_moon';

a_moon = mu_moon .* (r_sc_m / norm(r_sc_m)^3 - r_e_m / norm(r_e_m)^3);



a_pert = aj2 + a_moon;
end
