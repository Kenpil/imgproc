clear all;
clc;

img = imread('brainAct1.jpg');
height = length(img(:,1,1));
width = length(img(1,:,1));
redThre = 100;
greenThre = 50;
blueThre = 50;
actN = 0;
for i = 1:height
    for j = 1:width
        if(img(i,j,1) > redThre && img(i,j,2) < greenThre && img(i,j,3) < blueThre)
            actN = actN + 1;
        end
    end
end
actPlaceM = zeros(actN,2);
tmpN = 1;
for i = 1:height
    for j = 1:width
        if(img(i,j,1) > redThre && img(i,j,2) < greenThre && img(i,j,3) < blueThre)
            actPlaceM(tmpN,1) = i;
            actPlaceM(tmpN,2) = j;
            tmpN = tmpN + 1;
        end
    end
end
resultImg = img;
for i = 1:actN
    resultImg(actPlaceM(i,1),actPlaceM(i,2),1) = 0;
    resultImg(actPlaceM(i,1),actPlaceM(i,2),2) = 0;
    resultImg(actPlaceM(i,1),actPlaceM(i,2),3) = 255;
end
imshow(resultImg);

%%

idx = kmeans(actPlaceM, 2);

for j = 1:2
    colorRand = randi([0 255],3,1);
    for i = 1:actN
        if(idx(i,1) == j)
            resultImg(actPlaceM(i,1),actPlaceM(i,2),1) = colorRand(1,1);
            resultImg(actPlaceM(i,1),actPlaceM(i,2),2) = colorRand(2,1);
            resultImg(actPlaceM(i,1),actPlaceM(i,2),3) = colorRand(3,1);
        end
    end
end
imshow(resultImg);

%%
BIC = bic(actPlaceM,2,idx)
% resultImg(89,196,2) = 255; 
% resultImg(111,291,2) = 255; 
% resultImg(221,244,2) = 255; 
% resultImg(207,96,2) = 255;
% resultImg(138,77,2) = 255; 
% resultImg(161,204,2) = 255; 
% resultImg(50,128,2) = 255; 
% imshow(resultImg);

%%
function [xK, xIndex] = xmeans(x)
    R = length(x(:)); 
    index = ones(R,1);
    bic = BIC(x,1,index);
end

function BIC = bic(x,K,index)
    R = length(x(:,1)); % 要素の数:x(要素数,次元)
    M = length(x(1,:)); % 要素の次元
    Ri = zeros(K,1);
    vali = zeros(K,1); % 分散
    mui = zeros(K,M); % i番目クラスタの平均座標
    
    for i = 1:R
        for j = 1:K
            if(index(i) == j)
                Ri(j,1) = Ri(j,1) + 1;
                mui(j,:) = mui(j,:) + x(i,:);
            end
        end
    end
    for i = 1:K
        mui(i,:) = mui(i,:) / Ri(i,1);
    end
 
    for i = 1:R
        for j = 1:K
            if(index(i) == j)
                tmp = 0;
                for k = 1:M
                    tmp = tmp + power(x(i,k)-mui(j,k),2);
                end
                vali(j) = vali(j) + tmp;
            end
        end
    end
    vali = vali / (R-K);
    val = sum(vali);
    
    RiLogRi = 0;
    for i = 1:K
        RiLogRi = RiLogRi + Ri(i)*log(R(i));
    end
    IDn = -R/2*log(2*pi) - R*M/2*log(val) - (R-K)/2 + RiLogRi - R*log(R);
    pj = K*(M+1);
    BIC = IDn - pj/2*log(R);
end