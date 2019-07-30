%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  函数功能：标量非线性系统扩展Kalman滤波问题
%  状态函数：X(k+1)=0.5X(k)+2.5X(k)/(1+X(k)^2)+8cos(1.2k) +w(k)
%  观测方程：Z（k）=X(k)^2/20 +v(k)
%  详细原理介绍及中文注释请参考：
%  《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EKF_for_One_Div_UnLine_System
T=50;
Q=10;
R=1;
w=sqrt(Q)*randn(1,T);
v=sqrt(R)*randn(1,T);
x=zeros(1,T);
x(1)=0.1;
y=zeros(1,T);
y(1)=x(1)^2/20+v(1);
for k=2:T
    x(k)=0.5*x(k-1)+2.5*x(k-1)/(1+x(k-1)^2)+8*cos(1.2*k)+w(k-1);
    y(k)=x(k)^2/20+v(k);
end
Xekf=zeros(1,T);
Xekf(1)=x(1);
Yekf=zeros(1,T);
Yekf(1)=y(1);
P0=eye(1);
for k=2:T
    Xn=0.5*Xekf(k-1)+2.5*Xekf(k-1)/(1+Xekf(k-1)^2)+8*cos(1.2*k);
    Zn=Xn^2/20;
    F=0.5+2.5 *(1-Xn^2)/(1+Xn^2)^2;
    H=Xn/10;
    P=F*P0*F'+Q;
    K=P*H'*inv(H*P*H'+R);
    Xekf(k)=Xn+K*(y(k)-Zn);
    P0=(eye(1)-K*H)*P;
end
Xstd=zeros(1,T);
for k=1:T
    Xstd(k)=abs( Xekf(k)-x(k) );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on;box on;
plot(x,'-ko','MarkerFace','g');
plot(Xekf,'-ks','MarkerFace','b');
figure
hold on;box on;
plot(Xstd,'-ko','MarkerFace','g');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
