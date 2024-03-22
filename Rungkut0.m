clc;                                 
clear all;

%Data
data    = xlsread('DataSistem'); 
dataT   = .000000001*data';

%Nilai Parameter
Lambda  = 1000;
teta    = 15;
alfa1   = .000361925;
alfa2   = .00589634;
alfa3   = .000015081;
gama    = .000204;
k       = .001241771;
delta   = .025;
xi      = .00005;
beta    = -0.1621*10^-6;                

h = 0.1;  tfinal = 10^6;
N = ceil(tfinal/h) ;   % Calculates upto ceil(tfinal/h)

f1 = @(t, U,E,V,C)  Lambda-k*U*V-alfa1*U+gama*E-beta*U*C+xi*C-teta*U;     
f2 = @(t, U,E,V,C)  k*U*V-(alfa2+gama)*E;
f3 = @(t, U,E,V,C)  (alfa2+gama)*E-delta*V;
f4 = @(t, U,E,V,C)  beta*U*C-alfa3*C-xi*C+teta*U;

%Inisialisasi
t(1)=0;
%U(1)=dataT(1,1);  E(1)=dataT(2,1); V(1)=dataT(3,1); C(1)=dataT(4,1);
U(1)=data(1,1);  E(1)=data(1,2); V(1)=data(1,3); C(1)=data(1,4);
   
%Iterasi
for i=1:N;  %update loop
    t(i+1) = t(i) + h; %update time

 %Tahap 1%
 K1U = f1( t(i) , U(i), E(i), V(i), C(i) );
 K1E = f2( t(i) , U(i), E(i), V(i), C(i) );
 K1V = f3( t(i) , U(i), E(i), V(i), C(i) );
 K1C = f4( t(i) , U(i), E(i), V(i), C(i) );

 %Tahap 2%
 K2U = f1( t(i) + 0.5*h,   U(i) + 0.5*h*K1U,   E(i) +  0.5*h*K1E,   V(i) + 0.5*h* K1V ,   C(i) + 0.5*h* K1C );
 K2E = f2( t(i) + 0.5*h,   U(i) + 0.5*h*K1U,   E(i) +  0.5*h*K1E,   V(i) + 0.5*h* K1V ,   C(i) + 0.5*h* K1C );
 K2V = f3( t(i) + 0.5*h,   U(i) + 0.5*h*K1U,   E(i) +  0.5*h*K1E,   V(i) + 0.5*h* K1V ,   C(i) + 0.5*h* K1C );
 K2C = f4( t(i) + 0.5*h,   U(i) + 0.5*h*K1U,   E(i) +  0.5*h*K1E,   V(i) + 0.5*h* K1V ,   C(i) + 0.5*h* K1C );

 %Tahap 3%
 K3U = f1( t(i) + 0.5*h,   U(i) + 0.5*h*K2U,   E(i) +  0.5*h*K2E,   V(i) + 0.5*h* K2V ,   C(i) + 0.5*h* K2C );
 K3E = f2( t(i) + 0.5*h,   U(i) + 0.5*h*K2U,   E(i) +  0.5*h*K2E,   V(i) + 0.5*h* K2V ,   C(i) + 0.5*h* K2C );
 K3V = f3( t(i) + 0.5*h,   U(i) + 0.5*h*K2U,   E(i) +  0.5*h*K2E,   V(i) + 0.5*h* K2V ,   C(i) + 0.5*h* K2C );
 K3C = f4( t(i) + 0.5*h,   U(i) + 0.5*h*K2U,   E(i) +  0.5*h*K2E,   V(i) + 0.5*h* K2V ,   C(i) + 0.5*h* K2C );
 
 %Tahap 4%
 K4U = f1( t(i) + h,   U(i) + h*K3U,   E(i) +  h*K3E,   V(i) + h* K3V ,   C(i) +  h*K3C );
 K4E = f2( t(i) + h,   U(i) + h*K3U,   E(i) +  h*K3E,   V(i) + h* K3V ,   C(i) +  h*K3C );  
 K4V = f3( t(i) + h,   U(i) + h*K3U,   E(i) +  h*K3E,   V(i) + h* K3V ,   C(i) +  h*K3C ); 
 K4C = f4( t(i) + h,   U(i) + h*K3U,   E(i) +  h*K3E,   V(i) + h* K3V ,   C(i) +  h*K3C );   
     

 U(i+1) = U(i) + (1/6)*(K1U+2*K2U+2*K3U+K4U)*h;  
 E(i+1) = E(i) + (1/6)*(K1E+2*K2E+2*K3E+K4E)*h; 
 V(i+1) = V(i) + (1/6)*(K1V+2*K2V+2*K3V+K4V)*h;
 C(i+1) = C(i) + (1/6)*(K1C+2*K2C+2*K3C+K4C)*h;
    
end

%Menampilkan grafik
figure(3);
plot(t,U,t,E,t,V,t,C, 'Linewidth',2)
legend('Penganguran', 'Pekerja', 'Pekerjaan', 'Penjahat')
xlabel('Tahun')
ylabel('Jumlah')