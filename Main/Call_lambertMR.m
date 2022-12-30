%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       ORBITAL MECHANICS                                 %
%                    Academic year 2020/2021                              %
%                                                                         %
%                  Module 3: Orbital manoeuvres                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% DESCRIPTION: Sample script illustating the use of lambertMR
%
% Camilla Colombo, 11/11/2016
% Juan Luis Gonzalo, 17/11/2020: Additional comments
%

clc;
clearvars;
close all;

% Read Help of lambertMR for details. You can open file lambertMR.m in the
% editor, or use command 'help lambertMR'
%
% lambertMR(RI,RF,TOF,MU,orbitType,Nrev,Ncase,optionsLMR)
%
% As inputs you need only:
%  - RI: Vector containing the initial position in Cartesian coordinates [L]
%  - RF: Vector containing the final position vector in Cartesian coordinates [L]
%  - TOF: Transfer time, time of flight [T]
%  - MU: Planetary constant of the primary [L^3/T^2]
% For the other parameters set:
%  - orbitType = 0: Direct orbit (0: direct; 1: retrograde)
%  - Nrev = 0: Zero-revolution case
%  - Ncase = 0: Not used for the zero-revolution case
%  - optionsLMR = 0: No display


% muSun = 0.39860e6;      % Sun's gravitational parameter [km^3/s^2];
% ToF = 50*86400;         % Time in [s];
% r1 = [10000,0,0];       % Initial position vector [km]
% r2 = [20000,100,0];     % Final position vector [km]

%%EXCERCISE 1
mu_e =astroConstants(13);
r1 = [-21800,37900,0];
r2 = [27300,27700,0];
ToF = 15 * 3600 + 6 * 60 + 40;


[a_t,p_t,e_t,ERROR,v1,v2,T_par,theta] = lambertMR( r1, r2, ToF, mu_e, 0, 0, 0, 2);


figure(1)
DisegnaOrbita(a_t,e_t,w1,i1,omega1)
hold on;
grid on;
earth_sphere;
DisegnaPunto(p1,e1,w1,i1,omega1,theta1)
DisegnaPunto(p2,e2,w2,i2,omega2,theta2)
legend('Orbit','Earth','P1','P2');
%% EXCERCISE 2
clc;
clearvars;
close all;

mu = astroConstants(13);
%Orbit 1
a1 = 12500;
e1 = 0;
p1 = a1 * (1 - e1^2);
i1 = 0;
omega1 = 0;
w1 = 0;
theta1 = deg2rad(120);
[r1,v1]=kp2rv(a1,e1,i1,omega1,w1,theta1,mu);


%Orbit 2
a2 = 9500;
e2 = 0.3;
i2 = 0;
w2 = 0;
p2 = a2 * (1 - e2^2);
omega2 = 0;
theta2 = deg2rad(250);
ToF = 3300;
tspan_arc = linspace (0,ToF,1000);
[r2,v2]=kp2rv(a2,e2,i2,omega2,w2,theta2,mu);

%Lambert
[a_t,p_t,e_t,ERROR,v1t,v2t,T_par,thetaT] = lambertMR( r1, r2, ToF, mu, 0, 0, 0, 2);
T_trans = 2 * pi * sqrt(a_t^3/mu);
tspan_neg = linspace(ToF,T_trans,1000);

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
y0 = [r1';v1t'];
y0_neg = [r2';v2t'];
 [Tarc, Yarc ] = ode113( @(t,y) ode_2bp(t,y,mu), tspan_arc, y0, options );
 [Tneg, Yneg ] = ode113( @(t,y) ode_2bp(t,y,mu), tspan_neg, y0_neg, options );
 rt = Yarc(:,1:3);
 vt = Yarc(:,4:6);

Dv1 = norm(v1t - v1);
Dv2 = norm(v2 - v2t);
Dvtot = abs(Dv1 + Dv2);

%Plots
 DisegnaOrbita(a1,e1,w1,i1,omega1);         %Orbit 1
 hold on;
 earth_sphere;
 DisegnaOrbita(a2,e2,w2,i2,omega2);         %Orbit 2
 DisegnaPunto(p1,e1,w1,i1,omega1,theta1);   %Point 1
 DisegnaPunto(p2,e2,w2,i2,omega2,theta2);   %Point 2
 plot3(Yneg(:,1),Yneg(:,2),Yneg(:,3),'--','linewidth',2);   %Unused arc
 plot3(Yarc(:,1),Yarc(:,2),Yarc(:,3),'linewidth',2);     %Transfer Arc
 grid on;
 legend('Orbit 1','Earth','Orbit 2','P1','P2','Unused Arc','Transfer Arc');

 %% EXCERCISE 3
 clc;
 close all; clearvars;
 %Departure

  Departure_min = [2003,4,1,0,0,0];
  T_depar_min = date2mjd2000(Departure_min);   %T Departure Min
  Departure_max = [2003,8,1,0,0,0];
  T_depar_max = date2mjd2000(Departure_max);   % T Departure Max
  mu_sun = astroConstants(4);

 %Arrival

    arrive_min = [2003,9,1,0,0,0];
    T_arrive_min = date2mjd2000(arrive_min);
    arrive_max = [2004,3,1,0,0,0];
    T_arrive_max = date2mjd2000(arrive_max);
    
i = 1;
time_dep = linspace(T_depar_min,T_depar_max,800);
for t = time_dep

    [kep_E,ksun] = uplanet(t, 3);
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
 time_arr = linspace(T_arrive_min,T_arrive_max,800);
for t = time_arr

    [kep_Mars,ksun] = uplanet(t, 4);
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

Departure = linspace(T_depar_min,T_depar_max,800);
Arrive = linspace(T_arrive_min,T_arrive_max,800);
for i=1:length(v_dep)

    t1(i) = days2s(Departure(i));

    for k = 1:length(v_arr)


        t2(k) = days2s(Arrive(k));

        ToF(k,i) = t2(k) - t1(i);

      [a_t(k,i),p_t(k,i),e_t(k,i),ERROR,v_1t(k,:,i),v_2t(k,:,i),T_tr,theta_tr] = lambertMR(r_dep(i,:),r_arr(k,:),ToF(k,i),mu_sun,0,0,0,2);

       Dv1(k,i) = norm(v_1t(k,:,i) - v_dep(i,:));
       Dv2(k,i) = norm(v_arr(k,:) - v_2t(k,:,i));
       Dvtot(k,i) = abs(Dv1(k,i)) + abs(Dv2(k,i));
    end
end

for i = 1:length(Dvtot)

    for j = 1:length(Dvtot)

        if (Dvtot(i,j)>10)
            Dvtot(i,j) = NaN;
            
        end
    end
end

figure(1)
contour(Dvtot);
colorbar;

[Dvmin,c_min,r_min,at_min,et_min,ToF_min,Opt_Departure,Opt_Arrive] = FoundMin(Dvtot,ToF,a_t,Departure, Arrive);

figure(2)
DisegnaOrbita(a_E,e_E,w_E,i_E,Om_E)
hold on;
grid on;
axis equal;
DisegnaOrbita(a_M,e_M,w_M,i_M,Om_M)
DisegnaPunto(p_E,e_E,w_E,i_E,Om_E,theta_E(index_dep))
DisegnaPunto(p_M,e_M,w_M,i_M,Om_M,theta_M(index_arr(index_dep)))
[Dvmin,col_min,row_min,at_min,et_min,ToF_min,Opt_Departure,Opt_Arrive,r1_min,r2_min,v1_min,v2_min] = FoundMin(Dvtot,ToF,...
                                                                                                  a_t,e_t,r_dep,r_arr,v_arr,v_dep,Departure, Arrive);
PlotArcs(mu_sun,r1_min,v1_min,r2_min,v2_min,ToF_min,at_min);

%%
clc;close all;clearvars;
 Departure_min = [2003,4,1,0,0,0];

  Departure_max = [2003,8,1,0,0,0];

 %Arrival

    arrive_min = [2003,9,1,0,0,0];
   
    arrive_max = [2004,3,1,0,0,0];




[Dvtot,ToF] = DvTot(Departure_min,Departure_max,arrive_min,arrive_max,3,4,9);
[m,index_arr] = min(Dvtot);
[MinDv,index_dep] = min(m);
opt_dep = Departure(index_dep);
opt_arr = Arrive(index_arr(index_dep));

Departure_MinDv = mjd20002date(opt_dep);
Arrive_MInDv = mjd20002date(opt_arr);

figure(2)
% DisegnaOrbita(a_E,e_E,w_E,i_E,Om_E)
% hold on;
% grid on;
% axis equal;
% DisegnaOrbita(a_M,e_M,w_M,i_M,Om_M)
% DisegnaPunto(p_E,e_E,w_E,i_E,Om_E,theta_E(index_dep))
% DisegnaPunto(p_M,e_M,w_M,i_M,Om_M,theta_M(index_arr))




%%

