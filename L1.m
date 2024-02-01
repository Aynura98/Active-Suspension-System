clear

% Physical parameters
ms = 318.5;    % kg
mu = 35.5;     % kg
ks = 27000 ; % N/m
ku = 228000; % N/m

Kc = 938;

A = [ 0 1 0 0; [-ks 0 ks 0]/ms;0 0 0 1;[ks 0 -ks-ku 0]/mu];
B=[0; Kc/ms; 0; -Kc/mu];
E=[0; 0; 0; ku/mu];
H=[[-ks 0 ks 0]/ms;1 0 -1 0;0 0 1 0;0 0 0 0];
F=[Kc/ms;0;0;1];
L=[0;0;-1;0];

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

Ag = [A zeros(4,3) E*Cd;
      Bz*H Az Bz*L*Cd;
      zeros(1,4) zeros(1,3) Ad];
Bg = [B; Bz*F; 0];
Eg = [E*Dd;Bz*L*Dd;Bd];
Cg = [Dz*H Cz Dz*L*Cd];
Dg =Dz*F;
Fg = Dz*L*Dd;

nx=8;
nu=1;
nw=1;
nz=4;

gamma=sdpvar(1,1);
X=sdpvar(nx,nx);
Y=sdpvar(nu,nx);

lambda = 0.01;
mu = 1;

F1 = [(Ag*X+Bg*Y)+(Ag*X+Bg*Y)'+lambda*X Eg;
Eg' -mu*eye(1)]<=0;
F2 = [lambda*X zeros(nx,nu) (Cg*X+Dg*Y)';
zeros(nu,nx) (gamma-mu)*1 Fg';
(Cg*X+Dg*Y) Fg gamma*eye(4)]>=0;
F3 = X>=0;

%%%%DEFINIZIONE LMI%%%%
F_total=F1+F2+F3;

%%%%SETTAGGGIO DI YALMIP IN MODALITà 'SILENZIOSA'%%%%
opts=sdpsettings('verbose',0,'solver','');
%%%%RISOLUZIONE%%%%
solvesdp(F_total,gamma,opts);
Kc=double(Y)*inv(double(X));
a = eig(Ag+lambda/2*eye(8));