function PlotArcs(mu,r1,v1,r2,v2,ToF,a)

T_trans = 2 * pi * sqrt(a^3/mu);
tspan_arc = linspace(0:ToF,1000);
tspan_neg = linspace(ToF,T_trans,1000);

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
Y0 = [r1';v1'];
Y0_neg = [r2';v2'];
[Tarc, Yarc ] = ode113( @(t,y) ode_2bp(t,y,mu), tspan_arc, Y0, options );
[Tneg, Yneg ] = ode113( @(t,y) ode_2bp(t,y,mu), tspan_neg, Y0_neg, options );

 plot3(Yneg(:,1),Yneg(:,2),Yneg(:,3),'--','linewidth',2);   %Unused arc
 plot3(Yarc(:,1),Yarc(:,2),Yarc(:,3),'linewidth',2);     %Transfer Arc
 grid on;