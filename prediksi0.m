clc
clear
close all
data    = xlsread('DataSistem'); 
dataT   = .000000001*data';
h = 0.001; % time-step
t = 0:h:500; % simulation time
Nt = numel(t);
y = zeros(4,Nt);
Lambda= 5*10^6; teta= 15; alfa1= .000361925; alfa2= .00589634; alfa3= .000015081;
gama= .000204; k= .001241771; delta= .025; xi= .00005; beta= 2.4830665553128*10^-6;
y(: ,1)=[dataT(1,1); dataT(2,1); dataT(3,1); dataT(4,1)]; % initial values
f = @(t,y) [Lambda-k*y(1)*y(3)-alfa1*y(1)+gama*y(2)-beta*y(1)*y(4)+xi*y(4)-teta*y(1); k*y(1)*y(3)-(alfa2+gama)*y(2); (alfa2+gama)*y(2)-delta*y(3); beta*y(1)*y(4)-alfa3*y(4)-xi*y(4)+teta*y(1)]; % ode function     

for i=1:Nt-1
    k1=h * f(t(i)     , y(: , i));
    k2=h * f(t(i)+h/2 , y(: , i)+ k1/2);
    k3=h * f(t(i)+h/2 , y(: , i)+ k2/2);
    k4=h * f(t(i)+h   , y(: , i)+ k3);
    y(: , i+1) = y(: , i) + 1/6*(k1 + 2*k2 + 2*k3 + k4);
end
figure(3)
plot(t,y(1,:),'LineWidth',2)
hold on
plot(t,y(2,:),'LineWidth',2)
hold on
plot(t,y(3,:),'LineWidth',2)
hold on
plot(t,y(4,:),'LineWidth',2)
xlabel('Waktu')
ylabel('Jumlah')
grid on
legend({'Pengangguran','Pekerja','Pekerjaan','Penjahat'})