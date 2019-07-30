%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 功能说明：基于观测距离，扩展卡尔曼滤波完成对目标状态估计
%  详细原理介绍及中文注释请参考：
%  《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sys,x0,str,ts]=EKF(t,x,u,flag)
global Zdist;
global Xekf;  
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
        save('Xekf','Xekf');
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
global Xekf; 
Xekf=[x0];
global P;
P=zeros(2,2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sys=mdlUpdate(t,x,u,Q,R)
global Zdist; 
global Xekf;
global P;
Zdist=[Zdist,u]; 
%——————————————————————————————————————
x0=0;y0=0; 
Xold=Xekf(:,length(Xekf(1,:)));
Xpre=ffun(Xold);

Zpre=hfun(Xpre);

F=[1,0;0.1*cos(0.1*Xpre(1,1)),1];
H=[(Xpre(1,1)-x0)/Zpre,(Xpre(2,1)-y0)/Zpre];

Ppre=F*P*F'+Q;

K=Ppre*H'*inv(H*Ppre*H'+R);

Xnew=Xpre+K*(u-Zpre);

P=(eye(2)-K*H)*Ppre;

Xekf=[Xekf,Xnew];
sys=Xnew; 
function sys=mdlOutputs(t,x,u)
sys = x; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
