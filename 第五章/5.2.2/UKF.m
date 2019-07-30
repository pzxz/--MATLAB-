%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  无迹Kalman滤波在目标跟踪中的应用
%  详细原理介绍及中文注释请参考：
%  《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UKF
clc;clear;
T=1; 
N=60/T; 
X=zeros(4,N); 
X(:,1)=[-100,2,200,20]; 
Z=zeros(1,N); 
delta_w=1e-3; 
Q=delta_w*diag([0.5,1]) ;
G=[T^2/2,0;T,0;0,T^2/2;0,T];
R=5; 
F=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];
x0=200; 
y0=300;
Xstation=[x0,y0]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w=sqrtm(R)*randn(1,N);
for t=2:N
    X(:,t)=F*X(:,t-1)+G*sqrtm(Q)*randn(2,1);
end
for t=1:N
    Z(t)=Dist(X(:,t),Xstation)+w(t);
end
L=4;
alpha=1;
kalpha=0;
belta=2;
ramda=3-L;
for j=1:2*L+1
    Wm(j)=1/(2*(L+ramda));
    Wc(j)=1/(2*(L+ramda));
end
Wm(1)=ramda/(L+ramda);
Wc(1)=ramda/(L+ramda)+1-alpha^2+belta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xukf=zeros(4,N);
Xukf(:,1)=X(:,1);
P0=eye(4);
for t=2:N
    xestimate= Xukf(:,t-1);
    P=P0;
    cho=(chol(P*(L+ramda)))';
    for k=1:L
        xgamaP1(:,k)=xestimate+cho(:,k);
        xgamaP2(:,k)=xestimate-cho(:,k);
    end
    Xsigma=[xestimate,xgamaP1,xgamaP2];
    Xsigmapre=F*Xsigma;
    Xpred=zeros(4,1);
    for k=1:2*L+1
        Xpred=Xpred+Wm(k)*Xsigmapre(:,k);
    end
    Ppred=zeros(4,4);
    for k=1:2*L+1
        Ppred=Ppred+Wc(k)*(Xsigmapre(:,k)-Xpred)*(Xsigmapre(:,k)-Xpred)';
    end
    Ppred=Ppred+G*Q*G';
    chor=(chol((L+ramda)*Ppred))';
    for k=1:L
        XaugsigmaP1(:,k)=Xpred+chor(:,k);
        XaugsigmaP2(:,k)=Xpred-chor(:,k);
    end
    Xaugsigma=[Xpred XaugsigmaP1 XaugsigmaP2];
    for k=1:2*L+1
        Zsigmapre(1,k)=hfun(Xaugsigma(:,k),Xstation);
    end
    Zpred=0;
    for k=1:2*L+1
        Zpred=Zpred+Wm(k)*Zsigmapre(1,k);
    end
    Pzz=0;
    for k=1:2*L+1
        Pzz=Pzz+Wc(k)*(Zsigmapre(1,k)-Zpred)*(Zsigmapre(1,k)-Zpred)';
    end
    Pzz=Pzz+R;

    Pxz=zeros(4,1);
    for k=1:2*L+1
        Pxz=Pxz+Wc(k)*(Xaugsigma(:,k)-Xpred)*(Zsigmapre(1,k)-Zpred)';
    end
    K=Pxz*inv(Pzz);
    xestimate=Xpred+K*(Z(t)-Zpred);
    P=Ppred-K*Pzz*K';
    P0=P;
    Xukf(:,t)=xestimate;
end
for i=1:N
    Err_KalmanFilter(i)=Dist(X(:,i),Xukf(:,i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on;box on;
plot(X(1,:),X(3,:),'-k.');
plot(Xukf(1,:),Xukf(3,:),'-r+');
legend('真实轨迹','UKF轨迹')
figure
hold on; box on;
plot(Err_KalmanFilter,'-ks','MarkerFace','r')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=Dist(X1,X2)
if length(X2)<=2
    d=sqrt( (X1(1)-X2(1))^2 + (X1(3)-X2(2))^2 );
else
    d=sqrt( (X1(1)-X2(1))^2 + (X1(3)-X2(3))^2 );
end
function [y]=hfun(x,xx)
y=sqrt((x(1)-xx(1))^2+(x(3)-xx(2))^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
