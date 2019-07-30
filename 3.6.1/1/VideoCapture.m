%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 程序说明：这是一个视频捕获并录制的程序
%  详细原理介绍及中文注释请参考：
%  《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VideoCapture
aviobj = avifile('myAVI');
obj=videoinput('winvideo',1,'YUY2_320x240');
preview(obj);                
T=100;
k=0;
while (k<T)
    frame=getsnapshot(obj);     
    subplot(1,2,1)
    imshow(frame);
    frameRGB=ycbcr2rgb(frame);     
    subplot(1,2,2);
    imshow(frameRGB);
    drawnow;
    aviobj=addframe(aviobj,frameRGB);
    flushdata(obj)
    k=k+1
end
aviobj = close(aviobj); 
delete(obj); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
