function [a,p,e,i,omega,w,theta]= rv2kp(r,v,mu)
%Update for singularities of zero inclination :line 14
%Circular orbit : line 37

E=0.5*(norm(v))^2-mu/norm(r);
a=-mu/(2*E);
h=cross(r,v);
R=cross(v,h);
e=R/mu-r/norm(r);
p = a * (1 - norm(e)^2);
k=[0;0;1];
hcosi=h(3);
i=acos(hcosi/norm(h));
if(i == 0)
    Nx = 1;
    Ny = 0;
    Nz = 0;
    N = [Nx ;Ny ;Nz];
    omega = 0;
else
    Kh=cross(k,h);
    N=Kh/norm(Kh);
    Nx=N(1);         % cos(omega)
    Ny=N(2);
    omega=acos(Nx);
end
if(Ny<0)
    if(omega>pi)
        if(omega<2*pi)
            omega=2*pi-omega;
            
        end
    end
end
if (norm(e) == 0)
    w = 0;
else
    w=acos((dot(e,N)/norm(e)));
end
if(e(3)<0)
    if(w>pi)
        if(w<2*pi)
            w=2*pi-omega;
        end
    end
end
costheta=(dot(r,e)/(norm(r)*norm(e)));
theta=acos(costheta);
vr=dot(r,v);

if(vr<0)
   
     theta=2*pi-theta;
else
     
           

end

e=norm(e);


end
            


