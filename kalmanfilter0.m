clc
clear
close all
format long;

%% Data
data    = xlsread('DataSistem'); 
dataT   = .000000001*data';

%% Nilai Parameter
Lambda  = 5*10^3;
teta    = 15;
alfa1   = .000361925;
alfa2   = .00589634;
alfa3   = .000015081;
gama    = .000204;
k       = .001241771;
delta   = .025;
xi      = .00005;
beta    =  0;
dt      = .001;
n       = 11;

R      = 10^-6;
Q      = eye(5,5)*.5;
%Q(1,1) = .01;
%Q(2,2) = .01;
%Q(3,3) = .01;
Q(4,4) = .0005;
%Q(5,5) = .01;

P      = eye(5,5)*.5;
%P(1,1) = .001;
%P(2,2) = .001;
%P(3,3) = .001;
P(4,4) = .0005;
%P(5,5) = .001;

H      = [1 0 0 0 0; 
          0 1 0 0 0; 
          0 0 1 0 0; 
          0 0 0 1 0];

X_0            = [dataT(1,1); dataT(2,1); dataT(3,1); dataT(4,1); beta];
X(:,1)         = X_0;
X_hat(:,1)     = X_0;
X_hat_min(:,1) = X_0;

for i = 1:n-1

    %% Titik Setimbang
    U = delta/k;
    E = (Lambda*k^2^alfa2 - Lambda*k*X_hat(5,i)*delta + Lambda*k^2*xi - alfa2*delta*k*alfa2 + alfa1*delta^2*X_hat(5,i) - alfa1*delta*k*xi - teta*delta*k*alfa3) / (alfa2*k*(k*alfa2 - X_hat(5,i)*delta + k*xi));
    V = (alfa2 + gama)*(Lambda*k^2*alfa3 - Lambda*k*X_hat(5,i)*delta + Lambda*k^2*xi - alfa1*delta*k*alfa3 + alfa1*delta^2*X_hat(5,i) - alfa1*delta*k*xi - teta*delta*k*alfa3) / (delta*alfa2*k*(k*alfa3 - X_hat(5,i)*delta + k*xi));
    C = teta*delta / (k*alfa3 - X_hat(5,i)*delta + k*xi);

    %% Matriks Augmented State, Matriks G, dan Matriks H
    A = zeros(5,5);
    A(1,1) = (-k*V-alfa1-X_hat(5,i)*C-teta)*dt+1;
    A(1,2) = gama*dt;
    A(1,3) = -k*U*dt;
    A(1,4) = (-X_hat(5,i)*U+xi)*dt;
    A(1,5) = -U*C*dt;
    A(2,1) = k*V*dt;
    A(2,2) = (-alfa2-gama)*dt+1;
    A(2,3) = k*U*dt;
    A(3,2) = (alfa2+gama)*dt;
    A(3,3) = -delta*dt+1;
    A(4,1) = (X_hat(5,i)*C+teta)*dt;
    A(4,4) = (X_hat(5,i)*U-alfa3-xi)*dt+1;
    A(4,5) = U*C*dt;
    A(5,5) = 1;

        % True process
        X(:,i+1) = A*X(:,i) + Q*randn(5,1);
        
        % Observer model
        % Prediction based on system dynamics
        X_hat_min(:,i+1) = A*X_hat(:,i);
        P_min            = A*P*A' + Q;
        
        % Update based on measurement
        S = H*P_min*H'+ R;
        K = P_min*H'*inv(S);
        X_hat(:,i+1) = X_hat_min(:,i+1)+K*(dataT(:,i+1)-H*X_hat_min(:,i+1));
        P = (eye(5) - K*H)*P_min;
end

%Menampilkan nilai beta
disp('Nilai estimasi beta')
disp(X_hat(5,:))

%Menampilkan MSE
rmse1 = rmse(X_hat(1,:),dataT(1,:));
rmse2 = rmse(X_hat(2,:),dataT(2,:));
rmse3 = rmse(X_hat(3,:),dataT(3,:));
rmse4 = rmse(X_hat(4,:),dataT(4,:));

disp('RMSE Pengangguran :')
disp(rmse1)
disp('RMSE Pekerja :')
disp(rmse2)
disp('RMSE Pekerjaan :')
disp(rmse3)
disp('RMSE Penjahat :')
disp(rmse4)
%vektor.rmse = [rmse1; rmse2; rmse3; rmse4];
%disp(vektor.rmse)

figure (1);
subplot(2,2,1)
plot(dataT(1,:), '--b', 'LineWidth', 4)
hold on;
plot(X_hat(1,:), 'k')
legend('Actual','Estimation')
xlabel('Time')
ylabel('Unemployment')

subplot(2,2,2)
plot(dataT(2,:), '--b', 'LineWidth', 4)
hold on;
plot(X_hat(2,:), 'k')
legend('Actual','Estimation')
xlabel('Time')
ylabel('Employment')

subplot(2,2,3)
plot(dataT(3,:), '--b', 'LineWidth', 4)
hold on;
plot(X_hat(3,:), 'k')
legend('Actual','Estimation')
xlabel('Time')
ylabel('Vacancies')

subplot(2,2,4)
plot(dataT(4,:), '--b', 'LineWidth', 4)
hold on;
plot(X_hat(4,:), 'k')
legend('Actual','Estimation')
xlabel('Time')
ylabel('Crime')

figure (2);
plot(X_hat(5,:), 'k')
legend('Beta')
title('Interaksi Penjahat dan Pengangguran');
xlabel('Waktu')
ylabel('Tingkat Penambahan')



