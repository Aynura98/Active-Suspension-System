function K = H8_robust_f( Ag1,Ag2,Bg1,Bg2,Eg1,Eg2,Cg1,Cg2,Dg1,Dg2,Fg1,Fg2 )

    %% Dimensions
    nx=8;
    nu=1;
    nw=1;
    nz=4;
    
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
    
end