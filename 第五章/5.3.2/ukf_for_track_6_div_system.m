%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  功能说明： UKF在目标跟踪中的应用
%  参数说明： 1、状态6维，x方向的位置、速度、加速度；
%                y方向的位置、速度、加速度；
%             2、观测信息为距离和角度；
%  详细原理介绍及中文注释请参考：
%  《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ukf_for_track_6_div_system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=6;
t=0.5;
Q=[1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 0.01 0 0 0;
    0 0 0 0.01 0 0;
    0 0 0 0 0.0001 0;
    0 0 0 0 0 0.0001];
R = [100 0;
    0 0.001^2];
f=@(x)[x(1)+t*x(3)+0.5*t^2*x(5);x(2)+t*x(4)+0.5*t^2*x(6);...
    x(3)+t*x(5);x(4)+t*x(6);x(5);x(6)];
h=@(x)[sqrt(x(1)^2+x(2)^2);atan(x(2)/x(1))];
s=[1000;5000;10;50;2;-4];
x0=s+sqrtm(Q)*randn(n,1);
P0 =[100 0 0 0 0 0;
    0 100 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 0.1 0;
    0 0 0 0 0 0.1];
N=50;
Xukf = zeros(n,N);
X = zeros(n,N);
Z = zeros(2,N);
for i=1:N
    X(:,i)= f(s)+sqrtm(Q)*randn(6,1);
    s = X(:,i);
end
ux=x0;
for k=1:N
    Z(:,k)= h(X(:,k)) + sqrtm(R)*randn(2,1);
    [Xukf(:,k), P0] = ukf(f,ux,P0,h,Z(:,k),Q,R);
    ux=Xukf(:,k);
end
for k=1:N
    RMS(k)=sqrt( (X(1,k)-Xukf(1,k))^2+(X(2,k)-Xukf(2,k))^2 );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
t=1:N;
hold on;box on;
plot( X(1,t),X(2,t), 'k-')
plot(Z(1,t).*cos(Z(2,t)),Z(1,t).*sin(Z(2,t)),'-b.')
plot(Xukf(1,t),Xukf(2,t),'-r.')
legend('实际值','测量值','ukf估计值');
xlabel('x方向位置/米')
ylabel('y方向位置/米')
figure
box on;
plot(RMS,'-ko','MarkerFace','r')
xlabel('t/秒')
ylabel('偏差/米')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,P]=ukf(ffun,X,P,hfun,Z,Q,R)
L=numel(X);
m=numel(Z);
alpha=1e-2;
ki=0;
beta=2;
lambda=alpha^2*(L+ki)-L;
c=L+lambda;
Wm=[lambda/c 0.5/c+zeros(1,2*L)];
Wc=Wm;
Wc(1)=Wc(1)+(1-alpha^2+beta);
c=sqrt(c);
Xsigmaset=sigmas(X,P,c); 
[X1means,X1,P1,X2]=ut(ffun,Xsigmaset,Wm,Wc,L,Q);   
[Zpre,Z1,Pzz,Z2]=ut(hfun,X1,Wm,Wc,m,R);
Pxz=X2*diag(Wc)*Z2';
K=Pxz*inv(Pzz);
X=X1means+K*(Z-Zpre);
P=P1-K*Pxz';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xmeans,Xsigma_pre,P,Xdiv]=ut(fun,Xsigma,Wm,Wc,n,COV)
LL=size(Xsigma,2);
Xmeans=zeros(n,1);
Xsigma_pre=zeros(n,LL);
for k=1:LL
    Xsigma_pre(:,k)=fun(Xsigma(:,k));
    Xmeans=Xmeans+Wm(k)*Xsigma_pre(:,k);
end
Xdiv=Xsigma_pre-Xmeans(:,ones(1,LL));
P=Xdiv*diag(Wc)*Xdiv'+COV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xset=sigmas(X,P,c)
A = c*chol(P)';
Y = X(:,ones(1,numel(X)));
Xset = [X Y+A Y-A];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
