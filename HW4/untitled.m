figure(1);

fc1 = fcontour(@(x,y) sin(x/100*pi).*cos((y+40)/100*pi)-.5,'LevelList',[0]);
xlim([-200, 200])
ylim([-200, 200])


figure(2);
fcontour(@(x,y) x.*sin(y) - y.*cos(x)-1,'LevelList',[0]);

xlim([-200, 200])
ylim([-200, 200])

figure(3);

fc3 = fcontour(@(x,y) exp(-(x/40/3).^2-(y/40/3).^2) + exp(-(x/40+2).^2-(y/40+2).^2)-.9,'LevelList',[0]);
xlim([-200, 200])
ylim([-200, 200])

figure(4);

fsurf(@(x,y) -(sin(x/100*pi).*cos((y+40)/100*pi)-.5));
xlim([-200, 200])
ylim([-200, 200])
zlim([-1,1])

[x,y] = meshgrid(linspace(-200,200,1000),linspace(-200,200,1000));
figure(5);
(-9*(200-x)+3*(200+x))/400;
(-8*(200-y)+5*(200+y))/400;
surf(x,y, sin(sin(x*8/200)+cos(y*8/200))-cos(sin(x.*y*8/200*8/200)+cos(x*8/200))); shading interp

xlim([-200, 200])
ylim([-200, 200])
zlim([-10,10])

figure(6);
fsurf(@(x,y) -(exp(-(x/40/3).^2-(y/40/3).^2) + exp(-(x/40+2).^2-(y/40+2).^2)-.8));
xlim([-200, 200])
ylim([-200, 200])
zlim([-1,1])