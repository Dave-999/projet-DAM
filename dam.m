function [] = dam()

% Partie PRESSIONS

gamma=1.3;
D=79.5*10^-3;
R=47.75*10^-3;
L=199.5*10^-3;
Vc=474*10^-6;
tau=19.5;
Vmax=(tau/(tau-1))*Vc;
beta=L/R;
P0=10^5;
h=4*pi/800;

theta=linspace(-2*pi,2*pi,801);
P=zeros(1,801);
C1=P0*Vmax^gamma;
i=1;

while(-2*pi<=theta(i) && theta(i)<=-pi)
    P(i)=P0;
    a=i;
    i=i+1;
end
deltatheta=40*pi/180; 
thetad=-15*pi/180;

while(-pi<=theta(i) && theta(i)<=thetad);   
P(i)=C1/(((Vc/2)*(1-cos(theta(i))+beta-sqrt(beta^2-(sin(theta(i))*sin(theta(i)))))+Vc/(tau-1))^gamma);
b=i;
i=i+1;
end

while (theta(i)<=thetad+deltatheta && theta(i)>=thetad)
K1= f(theta(i-1),P(i-1));
K2= f(theta(i-1)+h,P(i-1)+h*K1);

P(i)= P(i-1)+(h/2)*(K1+K2);
c=i;
i=i+1;
end
C2=P(i-1)*((Vc/2)*(1-cos(theta(i-1))+beta-sqrt(beta^2-(sin(theta(i-1))^2)))+Vc/(tau-1))^gamma;

while(thetad+deltatheta<=theta(i) && theta(i)<=pi);   
P(i)=C2/((((Vc/2)*(1-cos(theta(i))+beta-sqrt(beta^2-(sin(theta(i))^2)))+Vc/(tau-1)))^gamma);
d=i;
i=i+1;
end
Pech=P(i-1)

k=2.5;
while(pi<=theta(i) && theta(i)<2*pi)
    P(i)=P0+(Pech-P0)*exp(-k*theta(i-200));
    e=i;
    i=i+1;
end

for i=1:801
P(i)=P(i)/10^5;
end

figure();
plot(theta(1:a),P(1:a),'blue','linewidth',2);
hold on;
plot(theta(a:b),P(a:b),'red','linewidth',2);
hold on;
plot(theta(b:c),P(b:c),'green','linewidth',2);
hold on;
plot(theta(c:d),P(c:d),'yellow','linewidth',2);
hold on;
plot(theta(d:e),P(d:e),'black','linewidth',2);
hold off;

title('Pression dans le cylindre en fonction de l''angle du vilebrequin')
legend('Admission','Compression','Explosion','Détente','Echappement');
xlabel('theta[rad]');
ylabel('p[bar]');

% Partie FORCES

for i=1:801
P(i)=P(i)*10^5; % conversion en Pascal
end

Wnormale = 2500/60; % 2500rpm
Welevee = 4000/60; % 4000rpm

Mpiston = 304*10^-3; % masses possibles : 304g et 372.5g
Mbielle = 591.5*10^-3; % sans les vis : 591.5g

D = 79.5*10^-3; % 79.5mm
R = 47.75*10^-3; % 47.75mm

FpiedNormale = zeros(1,801);
FpiedElevee = zeros(1,801);
FteteNormale = zeros(1,801);
FteteElevee = zeros(1,801);

for i=1:801
    FpiedNormale(i) = pi*(D^2)/4*P(i)-Mpiston*R*Wnormale^2*cos(theta(i));
    FpiedElevee(i) = pi*(D^2)/4*P(i)-Mpiston*R*Welevee^2*cos(theta(i));
    FteteNormale(i) = -pi*(D^2)/4*P(i)+(Mpiston+Mbielle)*R*Wnormale^2*cos(theta(i));
    FteteElevee(i) = -pi*(D^2)/4*P(i)+(Mpiston+Mbielle)*R*Welevee^2*cos(theta(i));
end

figure();
plot(theta,FpiedNormale,'blue','linewidth',2);
hold on;
plot(theta,FteteNormale,'red','linewidth',2);
hold off;

title('Efforts sur la bielle en fonction de l''angle du vilebrequin (\omega=2500rpm)');
legend('Efforts sur le pied de bielle','Efforts sur la tête de bielle');
xlabel('theta[rad]');
ylabel('F[N]');

figure();
plot(theta,FpiedElevee,'blue','linewidth',2);
hold on;
plot(theta,FteteElevee,'red','linewidth',2);
hold off;

title('Efforts sur la bielle en fonction de l''angle du vilebrequin (\omega=4000rpm)');
legend('Efforts sur le pied de bielle','Efforts sur la tête de bielle');
xlabel('theta[rad]');
ylabel('F[N]');

end

function dPdtheta = f(th,p)
gamma=1.3;
D=79.5*10^-3;
R=47.75*10^-3;
L=199.5*10^-3;
Vc=474*10^-6;
tau=19.5;
Vmax=(tau/(tau-1))*Vc;
beta=L/R;
P0=10^5;
Qtot=Vmax*1650*(P0/(8.3145*298.15))*(32*21/100+28*78/100+40/100) % transformer Vmax en m³
h=4*pi/400;

deltatheta=40*pi/180;
thetad=-15*pi/180;
V=(Vc/2)*(1-cos(th)+beta-sqrt(beta^2-(sin(th))^2))+Vc/(tau-1);
dQdtheta=pi*Qtot*sin((pi*th-pi*thetad)/deltatheta)/(2*deltatheta);
dVdtheta=(Vc/2)*sin(th)*(1+cos(th)/(sqrt(beta^(2)-sin(th)^(2))));
dPdtheta=-gamma*p*dVdtheta/V+(gamma-1)*dQdtheta/V;

end