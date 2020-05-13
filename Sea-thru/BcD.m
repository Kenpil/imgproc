clear all;
clc;

load('data/BcRemovedM.mat');
load('data/depthMap5.mat');
height = length(depthMap(:,1));
width = length(depthMap(1,:));

%%

load('data/acM5000.mat');
% acM(1500:1550,1500:1550,:)
% figure;
% imshow(acM);
%%
load('data/depthMap5.mat');
height = length(acM(:,1,1));
width = length(acM(1,:,1));

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

f = 2;
BcDSampleN = sum(round(rangeRankNs*0.001));
x = zeros(BcDSampleN,1);
BcDFull = zeros(BcDSampleN,3);
tmpN = 0;
for k = 1:10
    kRangeM = zeros(rangeRankNs(1,k),2); % [row, col];
    num = 1;
    for i = 1:height
        for j = 1:width           
            if depthMap(i,j) > k && depthMap(i,j) < k+1
                kRangeM(num,1) = i;
                kRangeM(num,2) = j;
                num = num + 1;
            end
        end
    end
    idx = randi([1 rangeRankNs(1,k)],1,round(rangeRankNs(1,k)*0.001));
    for i = 1:round(rangeRankNs(k)*0.001)
        colN = kRangeM(idx(i),1);
        rowN = kRangeM(idx(i),2);
        x(tmpN+i) = depthMap(colN,rowN);
        BcDFull(tmpN+i,:) = -log(f*acM(colN,rowN,:))/depthMap(colN,rowN);
%         BcDFull(tmpN+i,:) = f*acM(colN,rowN,:)*exp(depthMap(colN,rowN));
%         rgbIm(colN,rowN,1) = 1;
%         rgbIm(colN,rowN,2) = 0;
%         rgbIm(colN,rowN,3) = 0;
    end
    tmpN = tmpN + round(rangeRankNs(1,k)*0.001);
end

%%

BcoeffVals = zeros(3,4);
for i = 2:length(x)
    for j = 1:3
        if isinf(BcDFull(i,j))
            BcDFull(i,j) = BcDFull(i-1,j);
        end
    end
end
for i = 1:3
    [Bfit, ~] = fit(x, BcDFull(:,i), 'exp2', 'lower', [0 -Inf 0 -Inf], 'upper', [Inf 0 Inf 0]);
    coeffvals = coeffvalues(Bfit)
%     BcoeffVals(i,:) = coeffvals;
%     figure;
%     scatter(x,BcDFull(:,i));
%     hold on;
%     scatter(x,Bfit(x));
%     hold off;
end

%%

z = (0:0.001:20)';
expBcDzMapVals = zeros(length(z),3);
%BcoeffVals(1,:) = BcoeffVals(2,:);
%BcoeffVals = [3.39405778579470 -0.346198463967847 0.109868883514865 -8.76771675453485e-09; 3.39405778579470 -0.346198463967847 0.109868883514865 -8.76771675453485e-09; 3.39405778579470 -0.346198463967847 0.109868883514865 -8.76771675453485e-09];

% A =  BcoeffVals(i, 1);
% B =  BcoeffVals(i, 2);
% C =  BcoeffVals(i, 3);
% D =  BcoeffVals(i, 4);
% BcDz =  A*exp(B*z)+C*exp(D*z);
% figure;
% plot(z, BcDz);
    
for i = 1:3
    A =  BcoeffVals(i, 1);
    B =  BcoeffVals(i, 2);
    C =  BcoeffVals(i, 3);
    D =  BcoeffVals(i, 4);
    expBcDzMapVals(:,i) =  exp((A*exp(B*z)+C*exp(D*z)).*z);
%     figure;
%     plot(z, expBcDzMapVals(:,i));
%     figure;
%     plot(z, A*exp(B*z)+C*exp(D*z));
%     ylim([0 5]);
end

%%
 
% WcCoeff = [0.654, -0.4217; 0.7502, -0.2710; 0.654, -0.4217];
% Wc = WcCoeff(:,1)*exp(WcCoeff(:,2)*5);
% load('data/BcRemovedM.mat');
JcIm = zeros(length(depthMap(:,1,1)),length(depthMap(1,:,1)),3);
depthval = 1;
Wc = [1.5; 2; 2];
Wc = [1; 1; 1];
for j = 1:height
    for k = 1:width
        depthval = depthMap(j,k);
        for i = 1:3
            if depthval < 20
                JcIm(j,k,i) = BcRemovedIm(j,k,i)*expBcDzMapVals(round(1000*depthval),i)/Wc(i);
            end
        end
    end
end
%JcIm(1400:1450,1000:200,:)
figure;
imshow(JcIm);
title('JcIm');