clear all;
clc;

load('BcRemovedM_2.mat');
load('depthMap4910.mat');
load('acM2095000.mat');
height = length(depthMap(:,1));
width = length(depthMap(1,:));
%grayIm = rgb2gray(BcRemovedIm);
grayIm = rgb2gray(acM);

%%
%figure;
%imshow(BcRemovedIm);
%figure;
%imshow(acM*0.75);
%acM(:,:,2) = 0;
%acM(:,:,3) = 0;
%figure;
%imshow(acM);
depthlim = [0 4];
imagesc(depthMap,depthlim);
%%

bandStart = 1;
bandWidth = 7900;
rangeLength = round((3.2-0.9)*1000);
rangeRankNs = zeros(rangeLength,5);
rangeRankNs(:,5) = linspace(0.9,3.2,rangeLength);
rowTmp = 0;
colTmp = 0;
tmp = 0;
for i = 1:height
    for j = 1:1500
    %for j = bandStart:bandStart+bandWidth
        tmp = round(1000*(depthMap(i,j)-0.9));
        if(tmp > 0 && tmp < rangeLength)
            for k = 1:3
                rangeRankNs(tmp,k) = rangeRankNs(tmp,k) + acM(i,j,k);
            end
            rangeRankNs(tmp,4) = rangeRankNs(tmp,4) + 1;
        end
    end
end
for i = 1:height
    for j = 6500:7900
        tmp = round(1000*(depthMap(i,j)-0.9));
        if(tmp > 0 && tmp < rangeLength)
            for k = 1:3
                rangeRankNs(tmp,k) = rangeRankNs(tmp,k) + acM(i,j,k);
            end
            rangeRankNs(tmp,4) = rangeRankNs(tmp,4) + 1;
        end
    end
end
for i = 1:3
    rangeRankNs(:,i) = rangeRankNs(:,i) ./ rangeRankNs(:,4); % 平均をとる
end
%figure;
%imshow(acM);

%%
% figure;
% plot(rangeRankNs(:,1));
% figure;
% plot(rangeRankNs(:,2));
% figure;
% plot(rangeRankNs(:,3));
%%
f = 2;
BcDFull = zeros(rangeLength,3);
bias = [0;0.5;0.5];
for i = 1:3
    %BcDFull(:,i) = -log(f*rangeRankNs(:,i))./(rangeRankNs(:,5)-bias(i,1));
    BcDFull(:,i) = -log(f*rangeRankNs(:,i))./(rangeRankNs(:,5));
    figure;
    plot(rangeRankNs(:,5), BcDFull(:,i));
end
BcDFull(isnan(BcDFull))=1;

BcoeffVals = zeros(3,4);
for i = 1:3
    [Bfit, ~] = fit(rangeRankNs(:,5), BcDFull(:,i), 'exp2', 'lower', [0 -Inf 0 -Inf], 'upper', [Inf 0 Inf 0]);
    %[Bfit, ~] = fit(x, BcDFull(:,i), 'exp2');
    coeffvals = coeffvalues(Bfit);
    BcoeffVals(i,:) = coeffvals;
end

%%
z = (0:0.0001:15)';
BcD = zeros(length(z),3);
for i = 1:3
    A =  BcoeffVals(i, 1);
    B =  BcoeffVals(i, 2);
    C =  BcoeffVals(i, 3);
    D =  BcoeffVals(i, 4);
    BcD(:,i) =  (A*exp(B*z)+C*exp(D*z));
    %figure;
    %plot(z(1:3000),BcD(1:3000,i));
%      figure;
%      plot(z, expBcDzMapVals(:,i));
%    figure;
%    plot(z, A*exp(B*z)+C*exp(D*z));
    ylim([0 3.3]);
end

%%
zHat = zeros(rangeLength,3);
for i = 1:3
    zHat(:,i) = -log(2*rangeRankNs(:,i))./BcD(round(rangeRankNs(:,5)*1000),i);
end
%%
for i = 1:3
    figure;
    scatter(rangeRankNs(:,5),zHat(:,i));
    ylim([0 3.2]);
end

%%

expBcDzMapVals = zeros(length(z),3);

%BcoeffVals(3,:) = BcoeffVals(2,:);

for i = 1:3
    A =  BcoeffVals(i, 1);
    B =  BcoeffVals(i, 2);
    C =  BcoeffVals(i, 3);
    D =  BcoeffVals(i, 4);
    expBcDzMapVals(:,i) =  exp((A*exp(B*z)+C*exp(D*z)).*z);
%      figure;
%      plot(z, expBcDzMapVals(:,i));
    %figure;
    %plot(z, A*exp(B*z)+C*exp(D*z));
    %ylim([-5 5]);
end

%%
 
% WcCoeff = [0.654, -0.4217; 0.7502, -0.2710; 0.654, -0.4217];
% Wc = WcCoeff(:,1)*exp(WcCoeff(:,2)*5);
% load('data/BcRemovedM.mat');
Wc = [8; 8; 8];
JcIm = zeros(length(depthMap(:,1,1)),length(depthMap(1,:,1)),3);
depthval = 1;
bias = [0,1,0];
for j = 1:height
    for k = 1:width
        depthval = depthMap(j,k);
        for i = 1:3
            if depthval < 20
                JcIm(j,k,i) = BcRemovedIm(j,k,i)*bias(i);%*expBcDzMapVals(round(10000*depthval),i);
                %JcIm(j,k,i) = bias(i)*expBcDzMapVals(round(1000*depthval),i);
            end
        end
    end
end
for i = 1:3
    %JcIm(:,:,i) = JcIm(:,:,i)/Wc(i);
end
%JcIm(1400:1450,1000:200,:)
figure;
imshow(JcIm);
title('JcIm');