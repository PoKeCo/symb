pkg load symbolic

syms d a k n t N X x eqnas F 
h=sym('h', [5 2]);
k=sym('k', [5 2]);

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
    B=[-1/2;+1/2;zeros(M+1-2,1)];
    B=[0;0;0;0;h(3,1);h(3,2);h(4,1);h(4,2);h(5,1);h(5,2)];
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

    Av(t)=c(row_var,col_var);
    Hv(t)=Av*Xans;

    H0(t)=[t, Hv(t)(1)];
    H1(t)=diff(H0(t),t);
    H2(t)=diff(H1(t),t);
    HKappaN(t)=H1(t)(1)*H2(t)(2)-H2(t)(1)*H1(t)(2);
    HKappaD(t)=(H1(t)(1)^2+H1(t)(2)^2)^(3/2);
    HKappa0(t)=HKappaN/HKappaD;
    HKappa1(t)=diff(HKappa0,t);
    HKappa2(t)=diff(HKappa1,t);
    hans=solve( HKappa0(-1)==k(3,1),
                HKappa0(+1)==k(3,2),
                HKappa1(-1)==k(4,1),
                HKappa1(+1)==k(4,2),
                HKappa2(-1)==k(5,1),
                HKappa2(+1)==k(5,2),h(3:5,1:2));
% $$$                 figure(1);
% $$$                 for i=1:m+1
% $$$                     subplot(m+1,1,i);
% $$$                     ezplot(Hv(t)(i),[-1 +1]);
% $$$                 end
                
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Basic Function F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
if( 1 )
    %F0(t)=[cos(t), sin(t)];
    %F0(t)=[t, Hv(t)(1)];
    %Ctrl=[-1.0,-1;
    %       0.0,-1; 
    %       0.0,+1; 
    %      +1.0,+1]; 
    
    Ctrl=[ 0,       0;
           7,       0; 
           10, (10-7); 
           10,     10]/10; 
    F0(t)=[(1-t)^3, 3*t*(1-t)^2, 3*t^2*(1-t), t^3]*Ctrl;
    F1(t)=diff(F0,t);
    F2(t)=diff(F1,t);


    FKappaN(t)=F1(t)(1)*F2(t)(2)-F2(t)(1)*F1(t)(2);
    FKappaD(t)=(F1(t)(1)^2+F1(t)(2)^2)^(3/2);
    FKappa0(t)=FKappaN/FKappaD;
    FKappa1(t)=diff(FKappa0,t);
    FKappa2(t)=diff(FKappa1,t);

    hsubs=-[sym(0),sym(0);
            sym(0),sym(0);
            FKappa0(0),FKappa0(+1);
            FKappa1(0),FKappa1(+1);
            3*FKappa0(0)+FKappa2(0),  (3*FKappa0(+1)+FKappa2(+1))];
    
    Hth0=Av*Xans;
    Hth0s(t)=subs(Hth0,h, hsubs);
   
    figure(1);
    for i=1:m+1
        subplot(m+1,1,i);
        ezplot(Hth0s(t)(i),[-1 +1]);
    end

    FND1(t)=sqrt(F1(t)(1)^2 + F1(t)(2)^2);
    FN1(t)=F1/FND1;

    G0(t)=[F0(t)(1)-FN1(t)(2)*Hth0s(t*2-1)(1), F0(t)(2)+FN1(t)(1)*Hth0s(t*2-1)(1)];
    G1(t)=diff(G0,t);
    G2(t)=diff(G1,t);


    if( 1 )
        
        GKappaN(t)=G1(t)(1)*G2(t)(2)-G2(t)(1)*G1(t)(2);
        GKappaD(t)=(G1(t)(1)^2+G1(t)(2)^2)^(3/2);
        GKappa0(t)=GKappaN/GKappaD;
        
        figure(2);
        hold off;
        ezplot(F0(t)(1),F0(t)(2),[0 +1]);
        hold on;
        ezplot(G0(t)(1),G0(t)(2),[0 +1]);
        hold off;
        
        figure(3);
        hold off;
        ezplot(FKappa0(t),[0 +1]);
        hold on;
        ezplot(GKappa0(t),[0 +1]);
        hold off;

        figure(4);
        hold off;
        ezplot(FKappa1(t),[0 +1]);
        hold on;
        GKappa1(t)=diff(GKappa0,t);
        ezplot(GKappa1(t),[0 +1]);
        hold off;

        GKappa2(t)=diff(GKappa1,t);
        figure(5);
        hold off;
        ezplot(FKappa2(t),[0 +1]);
        hold on;
        ezplot(GKappa2(t),[0 +1]);
        hold off;
        

    end
end