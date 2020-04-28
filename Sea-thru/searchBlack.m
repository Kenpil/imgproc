clear all;
clc;

rgbImOri = imread('image/rgbIm5.jpg');
%rgbIm = imread('image/grayscaleSample.jpg');
load('data/depthMap5.mat');
width = 2737;
height = 1827;
%  width = 1000;
%  height = 1000;
grayImOri = rgb2gray(rgbImOri);
rgbIm = im2double(rgbImOri);
grayIm = im2double(grayImOri);
figure;
imshow(grayIm);
figure;
depthlim = [2 10];
imagesc(depthMap,depthlim);
imwrite(grayIm,'image/grayIm.jpg');

blackN = round(width * height * 0.01);
maxDepth = 50;
maxBright = 0.45;

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
% for i = 1:blackN
%     colN = darkCandidateM(1,I(i));
%     rowN = darkCandidateM(2,I(i));
%     rgbIm(colN,rowN,1) = 255;
%     rgbIm(colN,rowN,2) = 0;
%     rgbIm(colN,rowN,3) = 0;
% end   
figure;
imshow(rgbIm);
%imwrite(rgbIm,'image/rgbImRed.jpg');
%%

% rgbLens = zeros(1000,1000);
% for i = 1:1000
%     for j = 1:1000
%         rgbLens(i,j) = rgbIm(i, j);
%     end
% end

x = zeros(blackN,1);
for i = 1:blackN
    x(i,1) = depthMap(darkCandidateM(1,I(i)), darkCandidateM(2,I(i)));
end 
BdataFull = zeros(blackN,3);
for j = 1:3
    for i = 1:blackN
        BdataFull(i,j) = rgbIm(darkCandidateM(1,I(1,i)), darkCandidateM(2,I(1,i)), j);
    end
end

BcoeffVals = zeros(3,4);
bEstimate = fittype(@(A,B,C,D,x) A*(1-exp(-1*B*x))+C*exp(-1*D*x));
for i = 1:3
    [Bfit, ~] = fit(x, BdataFull(:,i), bEstimate, 'StartPoint', [0.5, 2.5, 0.5, 2.5]);
    coeffvals = coeffvalues(Bfit)
    BcoeffVals(i,:) = coeffvals;
    figure;
    scatter(x,BdataFull(:,i));
    hold on;
    scatter(x,Bfit(x));
    hold off;
end