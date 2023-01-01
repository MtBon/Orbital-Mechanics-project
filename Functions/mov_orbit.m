function mov_orbit(Kep,T)

a = Kep(:,1);
e = Kep(:,2);
incl = Kep(:,3);
Omega = Kep(:,4);
omega = Kep(:,5);
k = 1;

h = figure;
%h.Visible = 'off';

 
 axis equal;
 rotate3d;

for i = 1 :50 : length(Kep)
 day = s2days(T(i));
 DisegnaOrbita(a(i) , e(i), omega(i), incl(i), Omega(i));
 hold on;
 earth_sphere;
title(num2str(day,'%.0f') ,'Days');
 
 drawnow;
 
 pause(1/100000);
 clf;

 movievector(k) = getframe(h,[20 20 400 400]);
 k = k + 1;

end

% myWriter = VideoWriter('Orbit evolution','MPEG-4');
% myWriter.FrameRate = 20;
% open(myWriter);
% writeVideo(myWriter,movievector);
% close (myWriter);

movie(movievector);

