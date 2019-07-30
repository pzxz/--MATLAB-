%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 程序说明：练习使用plot函数，线型设置
% 详细原理介绍请参考：
% 《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function main  % 主函数
A1=randn(1,10);
A2=randn(1,10);
A3=randn(1,10);
% 画图1
figure 
box on
hold on; %在同一个figure中多次调用plot，需要hold
plot(A1,'-r') %红色的实线
plot(A2,'-.g') %绿色的点画线
plot(A3,'-b.') %蓝色的实现，数据点为点
xlabel('X-axis')
ylabel('Y-axis')
% 画图2
figure
box on
hold on; %在同一个figure中多次调用plot，需要hold
plot(A1,'-ko','MarkerFaceColor','r') %黑色实线，红色圆圈数据点
plot(A2,'-cd','MarkerFaceColor','g') %蓝绿色实线，绿色菱形数据点
plot(A3,'-bs','MarkerFaceColor','b') %蓝色实线，蓝色方形数据点
% 画图3
figure
box on
hold on; %在同一个figure中多次调用plot，需要hold
%黑色实线，红色圆圈数据点,原点的颜色设为红色，大小为10，线宽为10，后三者顺序随便
plot(A1,'-ko','MarkerFaceColor','r','MarkerSize',10,'LineWidth',10) 
%蓝绿色实线，菱形数据点,原点的颜色设为绿色，大小为10，线宽为5，后三者顺序随便
plot(A2,'-cd','MarkerFaceColor','g','LineWidth',5,'MarkerSize',10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
