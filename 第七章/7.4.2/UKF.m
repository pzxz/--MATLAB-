%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 功能说明：基于观测距离，无迹卡尔曼滤波完成对目标状态估计
%  详细原理介绍及中文注释请参考：
%  《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sys,x0,str,ts]=UKF(t,x,u,flag)
global Zdist; 
global Xukf;
Q=diag([0.01,0.04]); 
R=1;
switch flag
    case 0 
        [sys,x0,str,ts]=mdlInitializeSizes;
    case 2 
        sys=mdlUpdate(t,x,u,Q,R);
    case 3 
        sys=mdlOutputs(t,x,u);
    case {1,4}
        sys=[];
    case 9 
        save('Xukf','Xukf');
        save('Zdist','Zdist');
    otherwise 
        error(['Unhandled flag = ',num2str(flag)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sys,x0,str,ts]=mdlInitializeSizes(N)
sizes = simsizes;
sizes.NumContStates  = 0; 
sizes.NumDiscStates  = 2;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 1; 
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0  = [0,0]';  
str = []; 
ts  = [-1 0];
global Zdist; 
Zdist=[];
global Xukf; 
Xukf=[x0];
global P;
P=0.01*eye(2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sys=mdlUpdate(t,x,u,Q,R)
global Zdist;
global Xukf;
global P;
Zdist=[Zdist,u];
%——————————————————————————————————————
Xin=Xukf(:,length(Xukf(1,:)));
[Xnew,P]=GetUkfResult(Xin,u,P,Q,R)
Xukf=[Xukf,Xnew];
sys=Xnew; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sys=mdlOutputs(t,x,u)
sys = x;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xout,P]=GetUkfResult(Xin,Z,P,Q,R) 
L=2; 
alpha=0.01;
kalpha=0;
belta=2;
ramda=alpha^2*(L+kalpha)-L;
for j=1:2*L+1
    Wm(j)=1/(2*(L+ramda));
    Wc(j)=1/(2*(L+ramda));
end
Wm(1)=ramda/(L+ramda);
Wc(1)=ramda/(L+ramda)+1-alpha^2+belta;

xestimate=Xin;

P
cho=(chol(P*(L+ramda)))';
for k=1:L
    xgamaP1(:,k)=xestimate+cho(:,k);
    xgamaP2(:,k)=xestimate-cho(:,k);
end
Xsigma=[xestimate,xgamaP1,xgamaP2]; 

for k=1:2*L+1
    Xsigmapre(:,k)=ffun(Xsigma(:,k));
end

Xpred=zeros(2,1);  
for k=1:2*L+1
    Xpred=Xpred+Wm(k)*Xsigmapre(:,k);
end
Ppred=zeros(2,2);  
for k=1:2*L+1
    Ppred=Ppred+Wc(k)*(Xsigmapre(:,k)-Xpred)*(Xsigmapre(:,k)-Xpred)';
end
Ppred=Ppred+Q;

chor=(chol((L+ramda)*Ppred))';
for k=1:L
    XaugsigmaP1(:,k)=Xpred+chor(:,k);
    XaugsigmaP2(:,k)=Xpred-chor(:,k);
end
Xaugsigma=[Xpred XaugsigmaP1 XaugsigmaP2];

for k=1:2*L+1   
    Zsigmapre(1,k)=hfun(Xaugsigma(:,k));
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
Pxz=zeros(2,1);
for k=1:2*L+1
    Pxz=Pxz+Wc(k)*(Xaugsigma(:,k)-Xpred)*(Zsigmapre(1,k)-Zpred)';
end

K=Pxz*inv(Pzz);             

Xout=Xpred+K*(Z-Zpred);         
P=Ppred-K*Pzz*K';             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
