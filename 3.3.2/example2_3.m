%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 卡尔曼滤波用于自由落体运动目标跟踪问题
% 详细原理介绍及中文注释请参考：
% 《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function example2_3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=1000; 
Q=[0,0;0,0];
R=1; 
W=sqrt(Q)*randn(2,N);
V=sqrt(R)*randn(1,N);
A=[1,1;0,1];
B=[0.5;1];
U=-1;
H=[1,0];
X=zeros(2,N);
X(:,1)=[95;1];
P0=[10,0;0,1];
Z=zeros(1,N);
Z(1)=H*X(:,1);
Xkf=zeros(2,N);
Xkf(:,1)=X(:,1);
err_P=zeros(N,2);
err_P(1,1)=P0(1,1);
err_P(1,2)=P0(2,2);
I=eye(2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=2:N
    X(:,k)=A*X(:,k-1)+B*U+W(k);
    Z(k)=H*X(:,k)+V(k);
    X_pre=A*Xkf(:,k-1)+B*U; 
    P_pre=A*P0*A'+Q;
    Kg=P_pre*H'*inv(H*P_pre*H'+R);
    Xkf(:,k)=X_pre+Kg*(Z(k)-H*X_pre);
    P0=(I-Kg*H)*P_pre;
    err_P(k,1)=P0(1,1);
    err_P(k,2)=P0(2,2);
end
messure_err_x=zeros(1,N);
kalman_err_x=zeros(1,N);
kalman_err_v=zeros(1,N);
for k=1:N
    messure_err_x(k)=Z(k)-X(1,k);
    kalman_err_x(k)=Xkf(1,k)-X(1,k);
    kalman_err_v(k)=Xkf(2,k)-X(2,k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(V);
title('messure noise')
figure
hold on,box on;
plot(messure_err_x,'-r.');
plot(kalman_err_x,'-g.');
legend('测量误差','kalman估计误差')
figureplot(kalman_err_v);
title('速度误差')
figure
plot(err_P(:,1));
title('位移误差均方值')
figure
plot(err_P(:,1));
title('速度误差均方值')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

