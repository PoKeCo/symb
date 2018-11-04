pkg load symbolic

syms d a k n t N X x f eqnas F

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve Offset Function H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( 1 )
    c(d,n)=factorial(n)/factorial(n-d)*t^(n-d);
    ct(d,n,t)=factorial(n)/factorial(n-d)*t^(n-d);
    M=9;
    m=(M-1)/2+1;

    row_ary=floor([0:M]'/2)*ones(1,M+1);
    col_ary=ones(M+1,1)*[M:-1:0];
    t_ary=repmat([-1;+1],(M+1)/2,1)*ones(1,M+1);
    %t_ary=repmat([0;+1],(M+1)/2,1)*ones(1,M+1);

    At=ct(row_ary,col_ary,t_ary);
    X=sym('X',[M+1 1]);
    B=[-0.5;+0.5;zeros(M+1-2,1)];
    B=[0;0;0;0;0;0;0;0;+1;+1];
    Xsol=solve(At*X==B,X);
    Xans=[Xsol.X1_1;
          Xsol.X2_1;
          Xsol.X3_1;
          Xsol.X4_1;
          Xsol.X5_1;
          Xsol.X6_1;
          Xsol.X7_1;
          Xsol.X8_1;
          Xsol.X9_1;
          Xsol.X10_1];

    row_var=[0:M]'*ones(1,M+1);
    col_var=ones(M+1,1)*[M:-1:0];

    Av=c(row_var,col_var);
    Hv(t)=Av*Xans;

    figure(1);
    for i=1:m+1
        subplot(m+1,1,i);
        ezplot(Hv(t)(i),[-1 +1]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Basic Function F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%F0(t)=[cos(t), sin(t)];
%F0(t)=[t, Hv(t)(1)];
%Ctrl=[-1.0,-1;
%       0.0,-1; 
%       0.0,+1; 
%      +1.0,+1]; 

Ctrl=[ 0.0, 0.0;
       0.7, 0.0; 
       1.0, (1-0.7); 
       1.0, 1.0]; 
F0(t)=[(1-t)^3, 3*t*(1-t)^2, 3*t^2*(1-t), t^3]*Ctrl;
F1(t)=diff(F0,t);
F2(t)=diff(F1,t);


FND1(t)=sqrt(F1(t)(1)^2 + F1(t)(2)^2);
FN1(t)=F1/FND1;

G0(t)=[F0(t)(1)-FN1(t)(2)*Hv(t*2-1)(1)*100, F0(t)(2)+FN1(t)(1)*Hv(t*2-1)(1)*100];
G1(t)=diff(G0,t);
G2(t)=diff(G1,t);


if( 1 )
    FKappaN(t)=F1(t)(1)*F2(t)(2)-F2(t)(1)*F1(t)(2);
    FKappaD(t)=(F1(t)(1)^2+F1(t)(2)^2)^(3/2);
    FKappa(t)=FKappaN/FKappaD;
    
    
    
    GKappaN(t)=G1(t)(1)*G2(t)(2)-G2(t)(1)*G1(t)(2);
    GKappaD(t)=(G1(t)(1)^2+G1(t)(2)^2)^(3/2);
    GKappa(t)=GKappaN/GKappaD;
    
    
    figure(2);
    hold off;
    ezplot(FKappa(t),[0 +1]);
    hold on;
    ezplot(GKappa(t),[0 +1]);
    hold off;

    figure(3);
    hold off;
    ezplot(F0(t)(1),F0(t)(2),[0 +1]);
    hold on;
    ezplot(G0(t)(1),G0(t)(2),[0 +1]);
    hold off;
end