clear all
clc

% Physical parameters
ms1 = 200;    % kg
ms2 = 400;    % kg
mu = 35.5;     % kg
ks = 27000 ; % N/m
ku = 228000; % N/m
Kc = 938;


A1 = [ 0 1 0 0; [-ks 0 ks 0]/ms1;0 0 0 1;[ks 0 -ks-ku 0]/mu];
B1=[0; Kc/ms1; 0; -Kc/mu];
E1=[0; 0; 0; ku/mu];
H1=[[-ks 0 ks 0]/ms1;1 0 -1 0;0 0 1 0;0 0 0 0];
F1=[Kc/ms1;0;0;1];
L1=[0;0;-1;0];

A2 = [ 0 1 0 0; [-ks 0 ks 0]/ms2;0 0 0 1;[ks 0 -ks-ku 0]/mu];
B2=[0; Kc/ms2; 0; -Kc/mu];
E2=[0; 0; 0; ku/mu];
H2=[[-ks 0 ks 0]/ms2;1 0 -1 0;0 0 1 0;0 0 0 0];
F2=[Kc/ms2;0;0;1];
L2=[0;0;-1;0];

p=1;

A = p*A1 + (1-p)*A2;
B = p*B1 + (1-p)*B2;
E = p*E1 + (1-p)*E2;
H = p*H1 + (1-p)*H2;
F = p*F1 + (1-p)*F2;
L = p*L1 + (1-p)*L2;


Wa = 0.1*tf(2*pi*50,[1 2*pi*50]);
We = 0.001*tf(10,1);
Wt = tf(2*pi*20,[1 2*pi*20]);
Wu = 0.001*tf([1 200],[1 1000000]);

Wd = tf(20, [1 20]);
Wd_s=ss(Wd);
[Ad,Bd,Cd,Dd]=ssdata(Wd_s);

Wz = append(Wa,We,Wt,Wu);
Wz_s=ss(Wz);
[Az,Bz,Cz,Dz]=ssdata(Wz_s);

Ag1 = [A1 zeros(4,3) E1*Cd;
      Bz*H1 Az Bz*L1*Cd;
      zeros(1,4) zeros(1,3) Ad];
Ag2 = [A2 zeros(4,3) E2*Cd;
      Bz*H2 Az Bz*L2*Cd;
      zeros(1,4) zeros(1,3) Ad];
Bg1 = [B1; Bz*F1; 0];
Bg2 = [B2; Bz*F2; 0];
Eg1 = [E1*Dd;Bz*L1*Dd;Bd];
Eg2 = [E2*Dd;Bz*L2*Dd;Bd];
Cg1 = [Dz*H1 Cz Dz*L1*Cd];
Cg2 = [Dz*H2 Cz Dz*L2*Cd];
Dg1 =Dz*F1;
Dg2 =Dz*F2;
Fg1 = Dz*L1*Dd;
Fg2 = Dz*L2*Dd;

nx=8;
nu=1;
nw=1;
nz=4;
%% Variables
X = sdpvar(nx,nx);
Y = sdpvar(nu,nx);

gamma=sdpvar(1);

V1=([(Ag1*X+Bg1*Y)+(Ag1*X+Bg1*Y)'        Eg1            (Cg1*X+Dg1*Y)';
                Eg1'             -gamma*eye(nw)        Fg1';
             (Cg1*X+Dg1*Y)               Fg1      -gamma*eye(nz)]<=0);
V2=([(Ag2*X+Bg2*Y)+(Ag2*X+Bg2*Y)'        Eg2            (Cg2*X+Dg2*Y)';
                Eg2'             -gamma*eye(nw)        Fg2';
             (Cg2*X+Dg2*Y)               Fg2     -gamma*eye(nz)]<=0);
V3=(X>=0);
V4=(gamma>=0);

V_total=V1+V2+V3+V4;

solvesdp(V_total,gamma);

K=double(Y)*inv(double(X));

