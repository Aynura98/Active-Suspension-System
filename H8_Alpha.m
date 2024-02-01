clear all

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

gamma=sdpvar(1);
X=sdpvar(nx);
Y=sdpvar(nu,nx);% LMI constraints

alfa=5;
teta=pi/4;
r=5;
%% Variables
    
S1=([2*alfa*X+Ag*X+X*Ag'+Bg*Y+Y'*Bg']<=0);
S2=([sin(teta)*(Ag*X+Bg*Y+X*Ag'+Y'*Bg') cos(teta)*(-Ag*X-Bg*Y+X*Ag'+Y'*Bg');
    cos(teta)*(Ag*X+Bg*Y-X*Ag'-Y'*Bg') sin(teta)*(Ag*X+Bg*Y+X*Ag'+Y'*Bg')]<=0);
S3=([-r*X  Ag;
      Ag' -r*X]<=0);

V1=([(Ag*X+Bg*Y)+(Ag*X+Bg*Y)'        Eg            (Cg*X+Dg*Y)';
                Eg'             -gamma*eye(nw)        Fg';
             (Cg*X+Dg*Y)               Fg      -gamma*eye(nz)]<=0);
V2=(X>=0);
V3=(gamma>=0);

V_total=V1+V2+V3+S1+S2+S3;
opts = sdpsettings;
opts.solver = 'lmilab';
solvesdp(V_total,gamma,opts)
K=double(Y)*inv(double(X));
