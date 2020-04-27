clear all;
clc;

rgbIm = imread('image/rgbIm2.jpg');
%rgbIm = imread('image/grayscaleSample.jpg');
load('Nachsholim/image_set_01/depthMap.mat');
width = 2737;
height = 1827;
% width = 370;
% height = 292;
grayIm = rgb2gray(rgbIm);
figure;
imshow(grayIm);
%imwrite(grayIm,'image/grayIm.jpg');

blackN = round(width * height * 0.01);
maxDepth = 50;
maxBright = 90;

darkCandidateN = 0;
for i = 1:height
    for j = 1:width
        if grayIm(i,j) < maxBright && depthMap(i,j) < maxDepth
            darkCandidateN = darkCandidateN + 1;
        end
    end
end
darkCandidateM = zeros(3,darkCandidateN);
tmp = 1;
for i = 1:height
    for j = 1:width
        if grayIm(i,j) < maxBright && depthMap(i,j) < maxDepth
            darkCandidateM(1,tmp) = i;
            darkCandidateM(2,tmp) = j;
            darkCandidateM(3,tmp) = grayIm(i,j);
            tmp = tmp + 1;
        end
    end
end

[~,I] = mink(darkCandidateM(3,:), blackN);
for i = 1:blackN
    colN = darkCandidateM(1,I(i));
    rowN = darkCandidateM(2,I(i));
    rgbIm(colN,rowN,1) = 255;
    rgbIm(colN,rowN,2) = 0;
    rgbIm(colN,rowN,3) = 0;
end

% for colN = 1300:1320
%     for rowN = 2500:2520
%         rgbIm(colN,rowN,1) = 0;
%         rgbIm(colN,rowN,2) = 255;
%         rgbIm(colN,rowN,3) = 0;
%     end
% end
    
figure;
imshow(rgbIm);
%imwrite(rgbIm,'image/rgbImRed.jpg');