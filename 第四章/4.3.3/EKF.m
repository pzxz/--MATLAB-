%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  扩展Kalman滤波在目标跟踪中的应用
%  详细原理介绍及中文注释请参考：
%  《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EKF
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
for t=2:N
    X(:,t)=F*X(:,t-1)+G*sqrtm(Q)*randn(2,1);
end
for t=1:N
    Z(t)=Dist(X(:,t),Xstation)+sqrtm(R)*randn;
end
Xekf=zeros(4,N);
Xekf(:,1)=X(:,1);
P0=eye(4);
for i=2:N
    Xn=F*Xekf(:,i-1);
    P1=F*P0*F'+G*Q*G';
    dd=Dist(Xn,Xstation);
    H=[(Xn(1,1)-x0)/dd,0,(Xn(3,1)-y0)/dd,0];
    K=P1*H'*inv(H*P1*H'+R);
    Xekf(:,i)=Xn+K*(Z(:,i)-dd);
    P0=(eye(4)-K*H)*P1;
end
for i=1:N
    Err_KalmanFilter(i)=Dist(X(:,i),Xekf(:,i));
end
figure
hold on;box on;
plot(X(1,:),X(3,:),'-k.');
plot(Xekf(1,:),Xekf(3,:),'-r+');
legend('真实轨迹','EKF轨迹')
figure
hold on; box on;
plot(Err_KalmanFilter,'-ks','MarkerFace','r')
function d=Dist(X1,X2);
if length(X2)<=2
    d=sqrt( (X1(1)-X2(1))^2 + (X1(3)-X2(2))^2 );
else
    d=sqrt( (X1(1)-X2(1))^2 + (X1(3)-X2(3))^2 );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%