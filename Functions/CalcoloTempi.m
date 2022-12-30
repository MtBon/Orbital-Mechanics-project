function [Dt12] = CalcoloTempi(a,e,theta1,theta2)
mu=398600.433;
E1=2*atan(tan(theta1/2)*sqrt((1-e)/(1+e)));
E2=2*atan(tan(theta2/2)*sqrt((1-e)/(1+e)));
T=2*pi*sqrt(a^3/mu);
if(theta1>pi)
    Eespl=2*atan(tan((2*pi-theta1)/2)*sqrt((1-e)/(1+e)));
    t1=T-sqrt(a^3/mu)*(Eespl-e*sin(Eespl));
end
if(0<=theta1)
    if(theta1<=pi)
t1=sqrt(a^3/mu)*(E1-e*sin(E1));
    end
end
if(theta2>pi)
    Eespl2=2*atan(tan((2*pi-theta2)/2)*sqrt((1-e)/(1+e)));
    t2=T-sqrt(a^3/mu)*(Eespl2-e*sin(Eespl2));
end
if(0<=theta2)
    if(theta2<=pi)
t2=sqrt(a^3/mu)*(E2-e*sin(E2));
    end
end
if (theta2>theta1)
    Dt12=t2-t1; 
else
    Dt12=t2-t1+T;
end

% i pedici vanno cambiati nello script del laboratorio in base alla fase in
% cui ci si trova

end