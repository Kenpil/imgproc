clear all;
clc;

load('BcRemovedM_2.mat');
load('depthMap4910.mat');
load('acM2095000.mat');
height = length(depthMap(:,1));
width = length(depthMap(1,:));
%grayIm = rgb2gray(BcRemovedIm);
grayIm = rgb2gray(acM);
%figure;
%imshow(acM)
%figure
%imshow(grayIm);
%%

bandStart = 500;
bandWidth = 100;
rangeRankNs = zeros(1,10);
kRange = linspace(0.9,3.2,11);
rowTmp = 0;
colTmp = 0;
for k = 1:10
    randTmp = [rand rand rand];
    for i = 1:height
        for j = bandStart:bandStart+bandWidth
            if depthMap(i,j) > kRange(k) && depthMap(i,j) < kRange(k+1) % km~(k+1)mの範囲に何ピクセルがあるか数える
                rangeRankNs(1,k) = rangeRankNs(1,k) + 1;
                %acM(i,j,:) = randTmp;
            end
        end
    end
end
%figure;
%imshow(acM);
%%
%load('acM2095000.mat');
%rgbCheck = cast(acM,'uint8');

f = 2;
usedClass = 0.4; % 上位1%を計算に使う
usedNum = 0.01; % usedClass個のうちusedNum割を計算に使う
whiteN = sum(round(round(rangeRankNs*usedClass)*usedNum));
tmpN = 0;
x = zeros(whiteN,1);
BcDFull = zeros(whiteN,3);
Ectmp = zeros(whiteN,4);
for k = 1:10
    kRangeM = zeros(rangeRankNs(1,k),3); % [x, y, depthv];
    num = 1;
    for i = 1:height
        for j = bandStart:bandStart+bandWidth        
            if depthMap(i,j) > kRange(k) && depthMap(i,j) < kRange(k+1)
                kRangeM(num,1) = i;
                kRangeM(num,2) = j;
                kRangeM(num,3) = grayIm(i,j);
                num = num + 1;
            end
        end
    end
    % まずはそれぞれのdepthごとに上位usedClass個の白いピクセルを取り出す
    % そのインデックスIを作製
    [~,I] = maxk(kRangeM(:,3), round(rangeRankNs(1,k)*usedClass));
    
    kRangeM2 = zeros(round(rangeRankNs(1,k)*usedClass),3);
    rowTmp = 0;
    colTmp = 0;
    %rand1Tmp = rand;
    %rand2Tmp = rand;
    %rand3Tmp = rand;
    for i = 1:round(rangeRankNs(1,k)*usedClass)
        kRangeM2(i,:) = kRangeM(I(i),:);
        %x(tmpN+i) = depthMap(rowN,colN);
        %BcDFull(tmpN+i,:) = -log(f*acM(rowN,colN,:))/x(tmpN+i);
%         BdataFull(tmpN+i,:) = -log(f*acM(rowN,colN,:))/x(tmpN+i); %* (1-exp(-0.17*depthMap(colN,rowN)));% expは謎の補正項
        %rowTmp = kRangeM2(i,1);
        %colTmp = kRangeM2(i,2);
        %grayIm(rowTmp,colTmp) = 1; 
        %acM(rowTmp,colTmp,1) = rand1Tmp;
        %acM(rowTmp,colTmp,2) = rand2Tmp;
        %acM(rowTmp,colTmp,3) = rand3Tmp;
    end
    clear I
    clear kRangeM
    %上位1%のうち下位usedNum
    [~,I2] = mink(kRangeM2(:,3), round(round(rangeRankNs(1,k)*usedClass)*usedNum));
    for i = 1:round(round(rangeRankNs(1,k)*usedClass)*usedNum)
        rowN = kRangeM2(I2(i),1);
        colN = kRangeM2(I2(i),2);
        x(tmpN+i) = depthMap(rowN,colN);
        BcDFull(tmpN+i,:) = -log(f*acM(rowN,colN,:))/x(tmpN+i);
        Ectmp(tmpN+i,1:3) = f*acM(rowN,colN,:);
        Ectmp(tmpN+i,4) = x(tmpN+i);
%         BdataFull(tmpN+i,:) = -log(f*acM(rowN,colN,:))/x(tmpN+i); %* (1-exp(-0.17*depthMap(colN,rowN)));% expは謎の補正項       
    end
    tmpN = tmpN + round(round(rangeRankNs(1,k)*usedClass)*usedNum);
end
BcDFull(isinf(BcDFull))=1;
%figure;
%imshow(acM);
%figure;
%imshow(grayIm);
clear grayIm

%%

BcoeffVals = zeros(3,4);
% for i = 2:length(x)
%     for j = 1:3
%         if isinf(BcDFull(i,j))
%             %BcDFull(i,j) = BcDFull(i-1,j);
%         end
%     end
% end
for i = 1:3
    %[Bfit, ~] = fit(x, BcDFull(:,i), 'exp2', 'lower', [0 -Inf 0 -Inf], 'upper', [Inf 0 Inf 0]);
    [Bfit, ~] = fit(x, BcDFull(:,i), 'exp2');
    coeffvals = coeffvalues(Bfit);
    BcoeffVals(i,:) = coeffvals;
    %figure;
    %scatter(x,BcDFull(:,i),'.');
    %hold on;
    %scatter(x,Bfit(x));
    %hold off;
end
%%
z = (0:0.001:15)';
BcD = zeros(length(z),3);
for i = 1:3
    A =  BcoeffVals(i, 1);
    B =  BcoeffVals(i, 2);
    C =  BcoeffVals(i, 3);
    D =  BcoeffVals(i, 4);
    BcD(:,i) =  (A*exp(B*z)+C*exp(D*z));
%      figure;
%      plot(z, expBcDzMapVals(:,i));
%    figure;
%    plot(z, A*exp(B*z)+C*exp(D*z));
%    ylim([-5 5]);
end
zHat = zeros(length(Ectmp(:,1)),3);
for i = 1:3
    for j = 1:length(Ectmp(:,1))
        zHat(j,i) = -log(Ectmp(j,i))./BcD(round(Ectmp(j,4)*1000));
    end
    figure;
    scatter(Ectmp(:,4),zHat(:,i));
end


%%

z = (0:0.001:15)';
expBcDzMapVals = zeros(length(z),3);

for i = 1:3
    A =  BcoeffVals(i, 1);
    B =  BcoeffVals(i, 2);
    C =  BcoeffVals(i, 3);
    D =  BcoeffVals(i, 4);
    expBcDzMapVals(:,i) =  exp((A*exp(B*z)+C*exp(D*z)).*z);
%      figure;
%      plot(z, expBcDzMapVals(:,i));
    figure;
    plot(z, A*exp(B*z)+C*exp(D*z));
    ylim([-5 5]);
end

%%
 
% WcCoeff = [0.654, -0.4217; 0.7502, -0.2710; 0.654, -0.4217];
% Wc = WcCoeff(:,1)*exp(WcCoeff(:,2)*5);
% load('data/BcRemovedM.mat');
Wc = [8; 8; 8];
JcIm = zeros(length(depthMap(:,1,1)),length(depthMap(1,:,1)),3);
depthval = 1;
for j = 1:height
    for k = 1:width
        depthval = depthMap(j,k);
        for i = 1:3
            if depthval < 20
                JcIm(j,k,i) = BcRemovedIm(j,k,i)*expBcDzMapVals(round(1000*depthval),i);
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