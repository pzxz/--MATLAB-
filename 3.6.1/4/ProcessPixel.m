%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 功能描述：操作视频帧中的像素，在每一帧中打上标签
%  详细原理介绍及中文注释请参考：
%  《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ProcessPixel
mov=aviread('C:\Program Files\MATLAB71\work\video.avi')
totalFrame=size(mov,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imageLogo=imread('1.png');
imageSize= imresize(imageLogo,1);
[height width channel] = size(imageSize);
figure('Name','Processing Pixel')
for i=1:totalFrame
    frameData=mov(i).cdata;
    subplot(1,2,1);
    imshow(frameData)
    xlabel('The original video')
    for ii=1:height
        for jj=1:width
            for kk=1:channel
                frameData(ii,jj,kk)=imageLogo(ii,jj,kk);
            end
        end
    end
    subplot(1,2,2);
    imshow(frameData)
    xlabel('The processed video')
    pause(0.02);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

