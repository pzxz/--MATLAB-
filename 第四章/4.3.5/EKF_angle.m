%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  扩展Kalman滤波在纯方位目标跟踪中的应用实例
%  详细原理介绍及中文注释请参考：
%  《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EKF_angle
clc;clear;
T=1;
N=40/T;
X=zeros(4,N);
X(:,1)=[0,2,1400,-10];
Z=zeros(1,N); 
delta_w=1e-4;
Q=delta_w*diag([1,1]) ;
G=[T^2/2,0;T,0;0,T^2/2;0,T];
R=0.1*pi/180;
F=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];
x0=0;
y0=1000; 
Xstation=[x0;y0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w=sqrtm(R)*randn(1,N);
for t=2:N
    X(:,t)=F*X(:,t-1)+G*sqrtm(Q)*randn(2,1);
end
for t=1:N
    Z(t)=hfun(X(:,t),Xstation)+w(t);
end
Xekf=zeros(4,N);
Xekf(:,1)=X(:,1);
P0=eye(4);
for i=2:N
    Xn=F*Xekf(:,i-1);
    P1=F*P0*F'+G*Q*G';
    dd=hfun(Xn,Xstation);
    D=Dist(Xn,Xstation);
    H=[-(Xn(3,1)-y0)/D,0,(Xn(1,1)-x0)/D,0];
    K=P1*H'*inv(H*P1*H'+R);
    Xekf(:,i)=Xn+K*(Z(:,i)-dd);
    P0=(eye(4)-K*H)*P1;
end
for i=1:N
  Err_KalmanFilter(i)=sqrt(Dist(X(:,i),Xekf(:,i)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on;box on;
plot(X(1,:),X(3,:),'-k.');
plot(Xekf(1,:),Xekf(3,:),'-r+');
legend('真实轨迹','EKF轨迹')
figure
hold on; box on;
plot(Err_KalmanFilter,'-ks','MarkerFace','r')
figure 
hold on;box on;
plot(Z/pi*180,'-r.','MarkerFace','r');
plot(Z/pi*180+w/pi*180,'-ko','MarkerFace','g');
legend('真实角度','观测角度');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cita=hfun(X1,X0)
if X1(3,1)-X0(2,1)>=0
    if X1(1,1)-X0(1,1)>0
        cita=atan(abs( (X1(3,1)-X0(2,1))/(X1(1,1)-X0(1,1)) ));
    elseif X1(1,1)-X0(1,1)==0
        cita=pi/2;
    else
        cita=pi/2+atan(abs( (X1(3,1)-X0(2,1))/(X1(1,1)-X0(1,1)) ));
    end
else
    if X1(1,1)-X0(1,1)>0
        cita=3*pi/2+atan(abs( (X1(3,1)-X0(2,1))/(X1(1,1)-X0(1,1)) ));
    elseif X1(1,1)-X0(1,1)==0
        cita=3*pi/2;
    else
        cita=pi+atan(abs( (X1(3,1)-X0(2,1))/(X1(1,1)-X0(1,1)) ));
    end
end
function d=Dist(X1,X2);
if length(X2)<=2
    d=( (X1(1)-X2(1))^2 + (X1(3)-X2(2))^2 );
else
    d=( (X1(1)-X2(1))^2 + (X1(3)-X2(3))^2 );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
