clear all

% Physical parameters
ms = 318.5;    % kg
mu = 35.5;     % kg
ks = 27000 ; % N/m
ku = 228000; % N/m

K = 938;

x0=[0.05;0.01;0.09;0.01];%initial condition

A = [ 0 1 0 0; [-ks 0 ks 0]/ms;0 0 0 1;[ks 0 -ks-ku 0]/mu];
B = [0 0; K/ms 0; 0 0; -K/mu ku/mu];

C = [[-ks 0 ks 0]/ms;1 0 -1 0;0 0 1 0];
D = [K/ms 0; 0 0;0 -1];

qcar = ss(A,B,C,D);

figure 
initial(qcar,x0);