function dPdtheta = f(th,p)
gamma=1.3;
D=79.5;
R=47.75;
L=199.5;
Vc=474;
tau=19.5;
Vmax=(tau/(tau-1))*Vc;
beta=L/R;
P0=100000;
Qtot=Vmax*1650/1000*(P0*32/(8.3145*298.15));
h=4*pi/400;
deltatheta=40*pi/180;
thetad=-15*pi/180;

V=(Vc/2)*(1-cos(th)+beta-sqrt(beta^2-(sin(th)^2)))+Vc/(tau-1);
dQdtheta=pi*Qtot*sin((pi*th-pi*thetad)/deltatheta)/(2*deltatheta);
dVdtheta=(Vc/2)*sin(th)*(1+cos(th)/(sqrt(beta^(2)-sin(th)^(2))));
dPdtheta=-gamma*p*dVdtheta/V+(gamma-1)*dQdtheta/V;

end