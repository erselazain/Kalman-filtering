clc
clear
close all
format long;

%% Nilai Parameter
Lambda = 5000;
teta = 15;
alfa1 = 0.000361925;
alfa2 = 0.00589634;
alfa3 = 0.000015081;
gama = 0.000204;
k = 0.001241771;
delta = 0.025;
xi = 0.00005;
beta = 0; %estimasi

%% Titik Setimbang
U_ = delta/k;
E_ = (Lambda*k^2^alfa2 - Lambda*k*beta*delta + Lambda*k^2*xi - alfa2*delta*k*alfa2 + alfa1*delta^2*beta - alfa1*delta*k*xi - teta*delta*k*alfa3) / (alfa2*k*(k*alfa2 - beta*delta + k*xi));
V_ = (alfa2 + gama)*(Lambda*k^2*alfa3 - Lambda*k*beta*delta + Lambda*k^2*xi - alfa1*delta*k*alfa3 + alfa1*delta^2*beta - alfa1*delta*k*xi - teta*delta*k*alfa3) / (delta*alfa2*k*(k*alfa3 - beta*delta + k*xi));
C_ = teta*delta / (k*alfa3 - beta*delta + k*xi);

J(1,1) = -k*V_-alfa1-beta*C_-teta;
J(1,2) = gama;
J(1,3) = -k*U_;
J(1,4) = -beta*U_+xi;
J(2,1) = k*V_;
J(2,2) = -alfa2-gama;
J(2,3) = k*U_;
J(2,4) = 0;
J(3,1) = 0;
J(3,2) = alfa2+gama;
J(3,3) = -delta;
J(3,4) = 0;
J(4,1) = beta*C_+teta;
J(4,2) = 0;
J(4,3) = 0;
J(4,4) = beta*U_-alfa3-xi;

%% Matriks Jacobian
%disp('Matriks J : ')
J = [J(1,1) J(1,2) J(1,3) J(1,4); J(2,1) J(2,2) J(2,3) J(2,4); J(3,1) J(3,2) J(3,3) J(3,4); J(4,1) J(4,2) J(4,3) J(4,4)];

%% Persamaan Karakteristik
%disp('Persamaan karakteristik : ')
L = poly(J);
%disp(L);

%% Nilai Eigen
disp('Akar-akar persamaan karakteristik (nilai eigen) : ')
K = eig(J);
disp(K);

%%Variabel
disp('Berikut nilai masing-masing Variabel : ')
disp('Pengangguran'); disp(U_);
disp('Pekerja');disp(E_);
disp('Pekerjaan');disp(V_);
disp('Penjahat');disp(C_);