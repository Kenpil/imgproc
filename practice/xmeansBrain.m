clear all;
clc;

img = imread('brainAct4.jpg');
KTrue = 4;
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

idx = kmeans(actPlaceM, KTrue);

for j = 1:KTrue
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

clear all;
clc;

pointLength = 420;
actPlaceM = zeros(pointLength,2);
KTrue = 4;
mu = [15 5; 15 -5; -15 5; -15 -5];
for i =1:KTrue
    sigma = [1 0; 0 1];
    R = mvnrnd(mu(i,:),sigma,pointLength/KTrue);
    actPlaceM(pointLength/KTrue*(i-1)+1:pointLength/KTrue*i,:) = R;
end
% plot(actPlaceM(:,1),actPlaceM(:,2),'.');
idx = kmeans(actPlaceM, KTrue);
gscatter(actPlaceM(:,1),actPlaceM(:,2),idx)

%%
clc;
[xK,xIndex] = xmeans(actPlaceM);
xK

%%
clc;
BIC = bic(actPlaceM,3,idx)
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
    R = length(x(:,1));
    M = length(x(1,:));
    index = ones(R,1);
    xK = 1;
    xIndex = ones(R,1);
    bicOri = bic(x,1,index)
    idxChild = kmeans(x, 2);
%     figure;
%     gscatter(x(:,1),x(:,2),idxChild);
    bicChild = bic(x,2,idxChild)
    if(bicOri < bicChild)
        xK = xK + 1;
        Rchild = zeros(2,1);
        for i = 1:R
            Rchild(idxChild(i)) = Rchild(idxChild(i)) + 1;
        end
        xChild1 = zeros(Rchild(1),M);
        xChild2 = zeros(Rchild(2),M);
        tmp1 = 1;
        tmp2 = 3;
        for i = 1:R
            if(idxChild(i) == 1)
                xChild1(tmp1,:) = x(i,:);
                tmp1 = tmp1 + 1;
            else
                xChild2(tmp2,:) = x(i,:);
                tmp2 = tmp2 + 1;
            end
        end
        [xKChild1,xIndexChild1] = xmeans(xChild1);
        [xKChild2,xIndexChild2] = xmeans(xChild2);
        xK = xK + xKChild1 - 1;
        xK = xK + xKChild2 - 1;
    end
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
    for i = 1:K
        vali(i) = vali(i) / (Ri(i)-K);
    end
    
    IDn = 0;
    for i = 1:K
        IDn = IDn - Ri(i)*log(2*pi)/2 - Ri(i)*M*log(vali(i))/2 - (Ri(i)-K)/2 + Ri(i)*log(Ri(i)) - Ri(i)*log(R);
    end
    mui
    R
    Ri
    vali
    A = - Ri(1)*log(2*pi)/2;
    B = - Ri(1)*M*log(vali(1))/2;
    C = - (Ri(1)-K)/2;
    D = Ri(1)*log(Ri(1));
    E = - Ri(1)*log(R);
    A
    B
    C
    D
    E
  
    pj = K*(M+1)*2;
    F = - pj/2*log(R)
    BIC = IDn - pj/2*log(R);
end