clear all;
clc;

load('data/BcRemovedM.mat');
load('data/depthMap5.mat');

%%

p = 0.3;
maxDepth = 50;
maxDepthThre = 0.05;
height = length(BcRemovedIm(:,1,1));
width = length(BcRemovedIm(1,:,1));
acM = zeros(height,width,3);
acMtmp = zeros(height,width,3);
repetition = 10;
for i=1:repetition
    i
    for j = 2:height-1
        for k = 2:width-1
            Ne = 0;
            actmp = zeros(3,1);
            if depthMap(j-1,k) < maxDepth && abs(depthMap(j-1,k)-depthMap(j,k)) < maxDepthThre
                Ne = Ne + 1;
                for l = 1:3
                    actmp(l,1) = actmp(l,1) + acM(j-1,k,l);
                end
            end
            if depthMap(j+1,k) < maxDepth && abs(depthMap(j+1,k)-depthMap(j,k)) < maxDepthThre
                Ne = Ne + 1;
                for l = 1:3
                    actmp(l,1) = actmp(l,1) + acM(j+1,k,l);
                end
            end
            if depthMap(j,k-1) < maxDepth && abs(depthMap(j,k-1)-depthMap(j,k)) < maxDepthThre
                Ne = Ne + 1;
                for l = 1:3
                    actmp(l,1) = actmp(l,1) + acM(j,k-1,l);
                end
            end
            if depthMap(j,k+1) < maxDepth && abs(depthMap(j,k+1)-depthMap(j,k)) < maxDepthThre
                Ne = Ne + 1;
                for l = 1:3
                    actmp(l,1) = actmp(l,1) + acM(j,k+1,l);
                end
            end
            if Ne > 0
                for l = 1:3
                    acMtmp(j,k,l) =  p*BcRemovedIm(j,k,l) + (1-p)*actmp(l,1)/Ne;
                end
            end
        end
    end
    acM = acMtmp;
end
save('data/acM.mat', 'acM');

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

f = 2;
BcDSampleN = sum(round(rangeRankNs*0.001));
x = zeros(BcDSampleN,1);
BcDFull = zeros(BcDSampleN,3);
tmpN = 0;
for k = 1:10
    kRangeM = zeros(rangeRankNs(1,k),2); % [x, y, depthv];
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
% imshow(rgbIm);

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
    BcoeffVals(i,:) = coeffvals;
    figure;
    scatter(x,BcDFull(:,i));
    hold on;
    scatter(x,Bfit(x));
    hold off;
end

%%


z = (0:0.001:20)';
expBcDzMapVals = zeros(length(z),3);
for i = 1:3
    A =  BcoeffVals(i, 1);
    B =  BcoeffVals(i, 2);
    C =  BcoeffVals(i, 3);
    D =  BcoeffVals(i, 4);
    expBcDzMapVals(:,i) =  exp((A*exp(B*z)+C*exp(D*z)).*z);
    figure;
    plot(z, expBcDzMapVals(:,i));
end

depthval = 1;
for j = 1:height
    for k = 1:width
        depthval = depthMap(j,k);
        for i = 1:3
            if depthval < 20
                acMtmp(j,k,i) = acM(j,k,i)*expBcDzMapVals(round(1000*depthval),i);
            end
        end
    end
end
figure;
imshow(acMtmp);