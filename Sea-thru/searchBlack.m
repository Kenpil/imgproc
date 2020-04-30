clear all;
clc;

%rgbImOri = imread('image/rgbIm3.jpg');
% rgbImOri = imread('D3/Raw/4910.ARW');
% imfinfo('image/rgbIm3.jpg')
% if true
%     row = 5700;
%     col = 7200;
%     imfinfo('D3/Raw/4910.ARW')
%     fin=fopen('D3/Raw/4910.ARW','r');
%     I=fread(fin,row*col,'uint8=>uint8');
%     Z=reshape(I,row,col);
%     Z=Z';
%     k=imshow(Z)
% end
rgbImOri = imread('Nachsholim/image_set_11/RGT.tif');
%rgbImOri = imread('image/rawimage.png');
%rgbIm = imread('image/grayscaleSample.jpg');
load('data/depthMap5.mat');
width = 2737;
height = 1827;
% width = 1454;
% height = 963;
grayImOri = rgb2gray(rgbImOri);
rgbIm = im2double(rgbImOri);
grayIm = im2double(grayImOri);
clear rgbImOri;
clear grayImOri;
% figure;
% imshow(grayIm);
%figure;
%depthlim = [0 10];
%imagesc(depthMap,depthlim);
% imwrite(grayIm,'image/grayIm.jpg');

%%

rangeRankNs = zeros(1,10);
for i = 1:height
    for j = 1:width
        for k = 1:10
            if depthMap(i,j) > k && depthMap(i,j) < k+1
                rangeRankNs(1,k) = rangeRankNs(1,k) + 1;
            end
        end
    end
end

blackN = sum(round(rangeRankNs*0.01));
x = zeros(blackN,1);
BdataFull = zeros(blackN,3);
tmpN = 0;
for k = 1:10
    kRangeM = zeros(rangeRankNs(1,k),3); % [x, y, depthv];
    num = 1;
    for i = 1:height
        for j = 1:width           
            if depthMap(i,j) > k && depthMap(i,j) < k+1
                kRangeM(num,1) = i;
                kRangeM(num,2) = j;
                kRangeM(num,3) = grayIm(i,j);
                num = num + 1;
            end
        end
    end
    [~,I] = mink(kRangeM(:,3), round(rangeRankNs(1,k)*0.01));
    for i = 1:round(rangeRankNs(k)*0.01)
        colN = kRangeM(I(i),1);
        rowN = kRangeM(I(i),2);
        x(tmpN+i) = depthMap(colN,rowN);
        BdataFull(tmpN+i,:) = rgbIm(colN,rowN, :) * (1-exp(-0.17*depthMap(colN,rowN)));% expは謎の補正項
%         rgbIm(colN,rowN,1) = 1;
%         rgbIm(colN,rowN,2) = 0;
%         rgbIm(colN,rowN,3) = 0;
    end
    tmpN = tmpN + round(rangeRankNs(1,k)*0.01);
end
% imshow(rgbIm);

%%

BcoeffVals = zeros(3,4);
bEstimate = fittype(@(A,B,C,D,x) A*(1-exp(-1*B*x))+C*exp(-1*D*x));
for i = 1:3
    [Bfit, ~] = fit(x, BdataFull(:,i), bEstimate, 'StartPoint', [0.5, 2.5, 0.5, 2.5], 'lower', [0 0 0 0], 'upper', [1 5 1 5]);
    coeffvals = coeffvalues(Bfit)
    BcoeffVals(i,:) = coeffvals;
    figure;
    scatter(x,BdataFull(:,i));
    hold on;
    scatter(x,Bfit(x));
    hold off;
end

%%

z = 0:0.01:20;
BcMapVals = zeros(3, length(z));
for i = 1:3
    A =  BcoeffVals(i, 1);
    B =  BcoeffVals(i, 2);
    C =  BcoeffVals(i, 3);
    D =  BcoeffVals(i, 4);
    BcMapVals(i,:) =  A*(1-exp(-1*B*z))+C*exp(-1*D*z);
%     figure;
%     plot(z, BcMapVals(i,:));
end
% 
% for i = 1:3
%     A =  BcoeffVals(i, 1);
%     B =  BcoeffVals(i, 2);
%     C =  BcoeffVals(i, 3);
%     D =  BcoeffVals(i, 4);
%     BcMapVals(i,:) =  A*(1-exp(-1*B*z))+C*exp(-1*D*z);
%     figure;
%     plot(z(1:800), BcMapVals(i,1:800));
% end

%%

BcIm = ones(height, width, 3);
depthval = 1;
for j = 1:height
    for k = 1:width
        depthval = depthMap(j,k);
        if depthval < 20
            BcIm(j,k,:) = BcMapVals(:,round(100*depthval))';
        end
    end
end
figure;
imshow(BcIm);

%%

BcRemovedIm = rgbIm - BcIm;
for i = 1:height
    for j = 1:width
        if BcRemovedIm(i,j) < 0
            BcRemovedIm(i,j) = 0;
        end
    end
end
figure;
imshow(BcRemovedIm);
imwrite(BcRemovedIm,'image/BcRemovedIm.jpg');