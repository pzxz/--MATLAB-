%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 功能说明：S函数仿真系统的状态方程
%  详细原理介绍及中文注释请参考：
%  《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sys,x0,str,ts]=SimuStateFunction(t,x,u,flag)
global Xstate;
switch flag
    case 0  
        [sys,x0,str,ts]=mdlInitializeSizes;
    case 2  
        sys=mdlUpdate(t,x,u);
    case 3  
        sys=mdlOutputs(t,x,u);
    case {1,4}
        sys=[];
    case 9   
        save('Xstate','Xstate');
    otherwise  
        error(['Unhandled flag = ',num2str(flag)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 0;  
sizes.NumDiscStates  = 2; 
sizes.NumOutputs     = 2; 
sizes.NumInputs      = 2;   
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;  
sys = simsizes(sizes);
x0  = [0,0]';          
str = [];             
ts  = [-1 0];  
global Xstate;
Xstate=[];
Xstate=[Xstate,x0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sys=mdlUpdate(t,x,u)
Xnew=ffun(x)+u;
sys=Xnew;
global Xstate;
Xstate=[Xstate,Xnew];
function sys=mdlOutputs(t,x,u)
sys = x;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
