%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 石油地震勘测输入白噪声估值器算法仿真程序
%  详细原理介绍及中文注释请参考：
%  《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Oil_Explore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
T=300;
F=[1,0;0.3,-0.5];
L=[-1,2]';
H=[1 1];
R=0.1;
n=2;
Qg=49;
longa=0.3;
Q=longa*Qg;
randn('seed',13)
g=sqrt(Qg)*randn(1,T+10);
rand('state',1);
para=rand(1,T+10);
for t=1:T+10
    if para(t)<longa
        b(t)=1;
    else
        b(t)=0;
    end
    w(t)=b(t)*g(t);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v=sqrt(R)*randn(1,T+10);
X=zeros(2,T+10);
Z=zeros(1,T+10);
Z(1)=H*X(:,1)+v(1);
for t=2:T+10
    X(:,t)=F*X(:,t-1)+L*w(t-1);
    Z(t)=H*X(:,t)+v(t);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P0=eye(n);
Xe=zeros(n,T+10);
PP=[];
for t=1:T+8
    XX=F*X(:,t);
    P=F*P0*F'+L*Q*L';
    PP=[PP,P];
    K(:,t)=P*H'*inv(H*P*H'+R);
    e(:,t)=Z(t)-H*XX;
    Xe(:,t)=XX+K(:,t)*e(:,t);
    P0=(eye(n)-K(:,t)*H)*P;
end
N=3;
for t=1:T+8
    Persai(:,:,t)=F*(eye(n)-K(:,t)*H);
    Qe(:,:,t)=H*PP(:,2*(t-1)+1:2*t)*H'+R;
end
for t=1:T+5
    M(1,t)=Q*L'*H'*inv(Qe(:,:,t+1));
    M(2,t)=Q*L'*Persai(:,:,t+1)'*H'*inv(Qe(:,:,t+2));
    M(3,t)=Q*L'*Persai(:,:,t+2)'*Persai(:,:,t+1)'*H'*inv(Qe(:,:,t+3));
end
for t=1:T
    wjian(1,t)=M(1,t+1)*e(t+1);
    wjian(2,t)=wjian(1,t)+M(2,t+2)*e(t+2);
    wjian(3,t)=wjian(2,t)+M(3,t+3)*e(t+3);
end
for Num=1:N
    subplot(3,1,Num);
    t=1:T;
    plot(t,wjian(Num,t),'b.');
    for t=1:T
        hh=line( [t,t],[0,w(t)] );
        set(hh,'color','k');
    end
    xlabel(['w(t)和',num2str(Num),'步平滑器'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

