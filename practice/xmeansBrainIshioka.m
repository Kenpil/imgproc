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
    
    if(K == 1)
        x_mu = zeros(R,2);
        for i = 1:M
            x_mu = x(:,i) - mui(i,1);
        end
        Vi = x_mu'*x_mu/R;
        ViInv = inv(Vi);
        l1 = 0;
        for i = 1:R
            l1 = l1 + x_mu(i,:)'*ViInv*x_mu(i,:);
        end
        l1 = -l1/2;
        BIC = -2*(-R*M*log(2*pi)/2 - R*log(det(Vi))/2 + l1) + M*(M+3)/2;
    end
    
    if(K == 2)
        x1_mu = zeros(Ri(1),2);
        x2_mu = zeros(Ri(2),2);
        tmp1 = 1;
        tmp2 = 1;
        for i = 1:R
            if(index(i) == 1)
                x1_mu(tmp1,:) = x(i,:) - mui(1,:);
                tmp1 = tmp1 + 1;
            else
                x2_mu(tmp1,:) = x(i) - mui(2,:);
                tmp2 = tmp2 + 1;
            end
        end
        Vi1 = x1_mu'*x1_mu/Ri(1);
        Vi1Inv = inv(Vi1);
        Vi2 = x2_mu'*x2_mu/Ri(2);
        Vi2Inv = inv(Vi2);
        
        normMu1_mu2 = pow((mui(1,1)-mui(2,1))^2 + (mui(1,2)-mui(2,2))^2, 0.5);
        beta = pow(normMu1_mu2/(det(Vi1)+det(Vi2)), 0.5);
        alpha = 0.5/normcdf(beta);
        l1 = 0;
        for i = 1:R
            l1 = l1 + x_mu(i,:)'*ViInv*x_mu(i,:);
        end
        l1 = -l1/2;
        BIC = -2*(-R*M*log(2*pi)/2 - R*log(det(Vi))/2 + l1) + M*(M+3)/2;
    end
end