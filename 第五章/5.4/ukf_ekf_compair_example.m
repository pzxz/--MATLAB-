%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 程序说明：对比UKF与EKF在非线性系统中应用的算法性能
%  详细原理介绍及中文注释请参考：
%  《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ukf_ekf_compair_example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=50;
L=1;
Q=10;
R=1;
W=sqrtm(Q)*randn(L,N);
V=sqrt(R)*randn(1,N);
X=zeros(L,N);
X(:,1)=[0.1]';
Z=zeros(1,N);
Z(1)=X(:,1)^2/20+V(1);
Xukf=zeros(L,N);
Xukf(:,1)=X(:,1)+sqrtm(Q)*randn(L,1);
Pukf=eye(L);
Xekf=zeros(L,N);
Xekf(:,1)=X(:,1)+sqrtm(Q)*randn(L,1);
Pekf=eye(L);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=2:N
    X(:,k)=0.5*X(:,k-1)+2.5*X(:,k-1)/(1+X(:,k-1)^2)+8*cos(1.2*k)+W(k);
    Z(k)=X(:,k)^2/20+V(k);
    [Xekf(:,k),Pekf]=ekf(Xekf(:,k-1),Pekf,Z(k),Q,R,k);
    [Xukf(:,k),Pukf]=ukf(Xukf(:,k-1),Pukf,Z(k),Q,R,k);
end
err_ekf=zeros(1,N);
err_ukf=zeros(1,N);
for k=1:N
    err_ekf(k)=abs(Xekf(1,k)-X(1,k));
    err_ukf(k)=abs(Xukf(1,k)-X(1,k));
end
XX=X-W;
err_ave_ekf=sum(err_ekf)/N
err_ave_ukf=sum(err_ukf)/N
figure
hold on;box on;
plot(X,'-r*');
plot(Xekf,'-ko');
plot(Xukf,'-b+');
legend('真实状态','EKF估计','UKF估计')
xlabel('时间k/s')
ylabel('状态值')
figure
hold on;box on;
plot(err_ekf,'-ro');
plot(err_ukf,'-b+');
xlabel('时间k/s')
ylabel('偏差绝对值')
legend('EKF估计','UKF估计')
function [Xout,Pout]=ekf(Xin,P,Zin,Q,R,k)
Xpre=0.5*Xin+2.5*Xin/(1+Xin^2)+8*cos(1.2*k);
F=[0.5+(2.5*(1+Xpre^2)-2.5*Xpre*2*Xpre)/(1+Xpre^2)^2];
Ppre=F*P*F'+Q;
Zpre=Xpre^2/20;
H=[Xpre/10];
K=Ppre*H'/(H*Ppre*H'+R);
Xout=Xpre+K*(Zin-Zpre);
Pout=(eye(1)-K*H)*Ppre;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xout,Pout]=ukf(X,P0,Z,Q,R,k)
L=1;
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
xestimate= X;
P=P0;
cho=(chol(P*(L+ramda)))';
for j=1:L
    xgamaP1(:,j)=xestimate+cho(:,j);
    xgamaP2(:,j)=xestimate-cho(:,j);
end
Xsigma=[xestimate,xgamaP1,xgamaP2];
for j=1:2*L+1
    Xsigmapre(:,j)=0.5*Xsigma(:,j)+2.5*Xsigma(:,j)/(1+Xsigma(:,j)^2)+8*cos(1.2*k);
end
Xpred=0;
for j=1:2*L+1
    Xpred=Xpred+Wm(j)*Xsigmapre(:,j);
end
Ppred=0;
for j=1:2*L+1
    Ppred=Ppred+Wc(j)*(Xsigmapre(:,j)-Xpred)*(Xsigmapre(:,j)-Xpred)';
end
Ppred=Ppred+Q;
chor=(chol((L+ramda)*Ppred))';
for j=1:L
    XaugsigmaP1(:,j)=Xpred+chor(:,j);
    XaugsigmaP2(:,j)=Xpred-chor(:,j);
end
Xaugsigma=[Xpred XaugsigmaP1 XaugsigmaP2];
for j=1:2*L+1
    Zsigmapre(1,j)=Xaugsigma(:,j)^2/20;
end
Zpred=0;
for j=1:2*L+1
    Zpred=Zpred+Wm(j)*Zsigmapre(1,j);
end
Pzz=0;
for j=1:2*L+1
    Pzz=Pzz+Wc(j)*(Zsigmapre(1,j)-Zpred)*(Zsigmapre(1,j)-Zpred)';
end
Pzz=Pzz+R;

Pxz=0;
for j=1:2*L+1
    Pxz=Pxz+Wc(j)*(Xaugsigma(:,j)-Xpred)*(Zsigmapre(1,j)-Zpred)';
end
K=Pxz*inv(Pzz);
xestimate=Xpred+K*(Z-Zpred);
P=Ppred-K*Pzz*K';
Pout=P;
Xout=xestimate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




