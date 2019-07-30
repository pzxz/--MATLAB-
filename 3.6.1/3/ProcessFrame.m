%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 功能描述：读取AVI，并将AVI视频的每一帧转为bmp图片存储
%  详细原理介绍及中文注释请参考：
%  《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ProcessFrame
mov=aviread('C:\Program Files\MATLAB71\work\video.avi')
totalFrame=size(mov,2);
figure('Name','show the movie');
movie(mov);
for i=1:totalFrame
    frameData=mov(i).cdata;
    bmpName=strcat('C:\Program Files\MATLAB71\work\imageFrame\image',...
        int2str(i),'.bmp');
    imwrite(frameData,bmpName,'bmp');
    pause(0.02);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
