function [E]=Ecc_anom_solve(t,e,a,mu,t0)
    n=sqrt(mu/a^3);
    E0=(n.*t0+(e.*sin(n.*t0))/(1-sin(n.*t0+e)+sin(n*t0)));
    for i=1:length(t)
        for j=1:length(e)
            Fun=@(Ec) n.*t(i)+e(j).*sin(Ec)-Ec;
            [E(i,j)]=fsolve(Fun,E0,optimoptions('fsolve','Display','off'));
        end
    end