function [Dvtot,a_t,e_t,Departure, Arrive] = DvTot(Departure_min,Departure_max,arrive_min,arrive_max,planet1,planet2,Dv_max,N)
        %This function compute the Porkchop plot for given windows of
        %departure and arrival. You can put a maximum Dv in order to neglet
        %unuseful Dv for your mission.
        %
        %INPUT:     -departure and arrival must be in Date form;
        %           -Dv_max : max Dv you want to see;
        %           -N : number of elements to for linspace(default 500)
        %           -Number associated to planets for Uplanet function. 
        %                   1:   Mercury
        %                   2:   Venus
        %                   3:   Earth
        %                   4:   Mars
        %                   5:   Jupiter
        %                   6:   Saturn
        %                   7:   Uranus
        %                   8:   Neptune
        %                   9:   Pluto
        %                   10:  Sun
        %
        %OUTPUT:    -Dvtot[NxN] = Cost of all possible manoeuvres;
        %           -ToF  = Time of flights;
        %
        %CALLED FUNCTION: -uplanet;
        %                 -astroConstants;
        %                 -date2mjd2000;
        %                 -kp2rv;
        %                 -days2s;
        %                 -lambertMR;
        %
if nargin==7
    N = 500;
end




%Departure

  T_depar_min = date2mjd2000(Departure_min);   %T Departure Min
  T_depar_max = date2mjd2000(Departure_max);   % T Departure Max
  mu_sun = astroConstants(4);

 %Arrival
  
    T_arrive_min = date2mjd2000(arrive_min);
    T_arrive_max = date2mjd2000(arrive_max);
    
i = 1;
time_dep = linspace(T_depar_min,T_depar_max,N);
for t = time_dep

    [kep_E,~] = uplanet(t, planet1);
    a_E = kep_E(1);
    e_E = kep_E(2);
    i_E = kep_E(3);
    p_E = a_E * (1-e_E^2);
    Om_E = kep_E(4);
    w_E = kep_E(5);
    theta_E(i,1) = kep_E(6);
   
    [r_dep(i,:),v_dep(i,:)]=kp2rv(a_E,e_E,i_E,Om_E,w_E,theta_E(i,1),mu_sun);
    i = i+1;

end

 i = 1;
 time_arr = linspace(T_arrive_min,T_arrive_max,N);
for t = time_arr

    [kep_Mars,~] = uplanet(t, planet2);
    a_M = kep_Mars(1);
    e_M = kep_Mars(2);
    p_M = a_M * (1 - e_M^2);
    i_M = kep_Mars(3);
    Om_M = kep_Mars(4);
    w_M = kep_Mars(5);
    theta_M(i,1) = kep_Mars(6);

    [r_arr(i,:),v_arr(i,:)]=kp2rv(a_M,e_M,i_M,Om_M,w_M,theta_M(i,1),mu_sun);
    i = i+1;
end

Departure = linspace(T_depar_min,T_depar_max,N);
Arrive = linspace(T_arrive_min,T_arrive_max,N);

for i=1:length(v_dep)

    t1(i) = days2s(Departure(i));

    for k = 1:length(v_arr)


        t2(k) = days2s(Arrive(k));

        ToF(k,i) = t2(k) - t1(i);

      [a_t(k,i),p_t(k,i),e_t(k,i),ERROR,v_1t(k,:,i),v_2t(k,:,i),T_tr,~] = lambertMR(r_dep(i,:),r_arr(k,:),ToF(k,i),mu_sun,0,0,0,2);

       Dv1(k,i) = norm(v_1t(k,:,i) - v_dep(i,:));
       Dv2(k,i) = norm(v_arr(k,:) - v_2t(k,:,i));
       Dvtot(k,i) = abs(Dv1(k,i)) + abs(Dv2(k,i));
    end
end

for i = 1:length(Dvtot)

    for j = 1:length(Dvtot)

        if (Dvtot(i,j)>Dv_max)
            Dvtot(i,j) = NaN;
            
        end
    end
end

figure(1)
contour(Dvtot);
colorbar;