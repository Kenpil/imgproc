clear all;
clc;

img = imread('brainAct7.jpg');
KTrue = 7;
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
figure;
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
figure;
imshow(resultImg);

%%

clear all;
clc;

pointLength = 1800; % 乱数生成する点の個数
actPlaceM = zeros(pointLength,2);
KTrue = 8; % 正解のクラスタの数

% pointLength個の点をKTrue個のクラスタになるよう乱数で求める
for i = 1:KTrue
    mu = [30*rand 30*rand];
    varRand = -0.05 + rand*0.1;
    rand1 = rand+0.1;
    rand2 = rand+0.1;
    sigma = [rand1 varRand; varRand rand2];
    X = mvnrnd(mu,sigma,pointLength/KTrue);
    actPlaceM((i-1)*pointLength/KTrue+1:i*pointLength/KTrue,:) = X;
end

% 乱数生成してみた点を表示
figure;
subplot(1,3,1);
scatter(actPlaceM(:,1),actPlaceM(:,2), '.','red');
title(['Original Points (',num2str(KTrue),' clusters)']);

% 正しいクラスタ数(KTrue個)でクラスタ数を明示して，従来のk-meansクラスタリングをしてみる
idx = kmeans(actPlaceM, KTrue);
% KTrue個のクラスタに応じた点の分布を図時
%figure;
subplot(1,3,2);
gscatter(actPlaceM(:,1),actPlaceM(:,2),idx);
title(['K-means Points (',num2str(KTrue),' clusters)']);

% 点のxy情報をmatファイルに保存
% save('C:\Users\Dell\WinRobots\develop\hw\biology\言語学C/actPlaceM.mat', 'actPlaceM');

% x-meansにより，クラスタ数を与えず，自動でクラスタ数を判別し，さらに点を分ける
% 第1引数はデータの（データ数,次元）行列
% xmeans()の呼び出し時には第2引数を1で初期化
% xkはクラスタリング結果のクラスタ数，xIndexはデータの属するクラスタ番号
[xK,xIndex] = xmeans(actPlaceM,1);
% 求めたクラスタ数
xK
% 求めたクラスタ数に応じた点の分布を図示
%figure;
subplot(1,3,3);
gscatter(actPlaceM(:,1),actPlaceM(:,2),xIndex);
title(['X-means Points (Auto Determined ',num2str(xK),' clusters)']);

%%

% xは(データ数, 次元)のデータ行列例えば(x,y)2次元分布なら(点1のx座標1,点1のy座標1; 点2のx座標1,点2のy座標1; ...; 点nのx座標1,点nのy座標1)となる．
% prepKは呼び出し時に1として初期化する．
% 返り値 xK: x-meansクラスタリング結果のクラスタ数
% 返り値 xIndex: x-meansクラスタリング結果の各データ点の属するクラスタ番号インデックス
function [xK, xIndex] = xmeans(x,preK)
    R = length(x(:,1)); % クラスタiの点数
    M = length(x(1,:)); % 点の次元，2次元平面上に点が散らばるならM=2
    index = ones(R,1); % 点がどのクラスタに属するか
    xK = 1;
    xIndex = ones(R,1);
    bicOri = bic(x,1,index);
    idxChild = kmeans(x, 2);
 
    bicChild = bic(x,2,idxChild);
    
    % bic(クラスタiそのまま) > bic(クラスタiを2分割した時) なら2分割した方が良い->再帰的に分割を続行
    if(bicOri > bicChild)
        Rchild = zeros(2,1);
        for i = 1:R
            Rchild(idxChild(i)) = Rchild(idxChild(i)) + 1;
        end
        xChild1 = zeros(Rchild(1),M);
        xChild2 = zeros(Rchild(2),M);
        tmp1 = 1;
        tmp2 = 1;
        for i = 1:R
            if(idxChild(i) == 1)
                xChild1(tmp1,:) = x(i,:);
                tmp1 = tmp1 + 1;
            else
                xChild2(tmp2,:) = x(i,:);
                tmp2 = tmp2 + 1;
            end
        end
        
        % 2分割後のクラスタそれぞれに対し，再帰的にx-meansを行い分割を続ける
        [xKChild1,xIndexChild1] = xmeans(xChild1,preK);
        [xKChild2,xIndexChild2] = xmeans(xChild2,preK + xKChild1);
        
        % 返り値を整理
        tmp1 = 1;
        tmp2 = 1;
        for i = 1:R
            if(idxChild(i) == 1)
                xIndex(i) = xIndexChild1(tmp1);
                tmp1 = tmp1 + 1;
            else
                xIndex(i) = xKChild1 + xIndexChild2(tmp2);
                tmp2 = tmp2 + 1;
            end
        end
        
        xK = xKChild1 + xKChild2;
        xIndex = xIndex;
    end
end

%% 


function BIC = bic(x,K,index)
    R = length(x(:,1)); % 要素の数:x(要素数,次元)
    M = length(x(1,:)); % 要素の次元
    Ri = zeros(K,1);
    mui = zeros(K,M); % i番目クラスタの平均座標，i番目のクラスタを2分割したら，(2,M)行列となり，[mu_x1 mu_y1; mu_x2 mu_y2]となる．
    BIC = 0;
    
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
    
    % クラスタi全体のBIC
    % BIC = -2log[クラスタiの最尤推定量] + M(M+3)log(クラスタiの点の数)/2
    if(K == 1)
        x_mu = zeros(R,2);
        % x - mu ベクトルx_mu
        for i = 1:M
            x_mu(:,i) = x(:,i) - mui(1,i);
        end
        Vi = x_mu'*x_mu/R;
        ViInv = inv(Vi);
        l1 = 0;
        % i全体の点のsum((x-mu)^(t)V^(-1)(x-mu))
        for i = 1:R
            l1 = l1 + x_mu(i,:)*ViInv*x_mu(i,:)';
        end
        l1 = -l1/2;
        BIC = -2*(-R*M*log(2*pi)/2 - R*log(det(Vi))/2 + l1) + M*(M+3)*log(R)/2;
    end
    
    
    % クラスタiを2分割した時のクラスタiのBIC
    % BIC = -2log[点ごとに2分割後に対応する平均・分散を用いたクラスタiの最尤推定量] + M(M+3)log(クラスタiの点の数)
    if(K == 2)
        x1_mu = zeros(Ri(1),2);
        x2_mu = zeros(Ri(2),2);
        tmp1 = 1;
        tmp2 = 1;
        % クラスタiの点がi_1,i_2のどちらに属すか判断し，x_mu = x - mu より平均から各点へのベクトルとして格納
        for i = 1:R
            if(index(i) == 1)
                x1_mu(tmp1,1) = x(i,1) - mui(1,1);
                x1_mu(tmp1,2) = x(i,2) - mui(1,2);
                tmp1 = tmp1 + 1;
            else
                x2_mu(tmp2,1) = x(i,1) - mui(2,1);
                x2_mu(tmp2,2) = x(i,2) - mui(2,2);
                tmp2 = tmp2 + 1;
            end
        end
        % i_1, i_2 それぞれの分散・共分散行列を求める
        Vi1 = x1_mu'*x1_mu/Ri(1);
        Vi1Inv = inv(Vi1);
        Vi2 = x2_mu'*x2_mu/Ri(2);
        Vi2Inv = inv(Vi2);
        
        normMu1_mu2_2 = power(mui(1,1)-mui(2,1),2) + power(mui(1,2)-mui(2,2),2);
        detVi1 = det(Vi1);
        detVi2 = det(Vi2);
        if detVi1 == 0 && detVi2 == 0 % 0割り対策
            beta = 0;
        else
            beta = power(normMu1_mu2_2/(detVi1+detVi2), 0.5);
        end
        alpha = 0.5/normcdf(beta);
        
        l1 = 0;
        l2 = 0;
        tmp1 = 1;
        tmp2 = 1;
        % i_1,i_2ごとのsum((x-mu)^(t)V^(-1)(x-mu))
        for i = 1:R
            if(index(i) == 1)
                l1 = l1 + x1_mu(tmp1,:)*Vi1Inv*x1_mu(tmp1,:)';
                tmp1 = tmp1 + 1;
            else
                l2 = l2 + x2_mu(tmp2,:)*Vi2Inv*x2_mu(tmp2,:)';
                tmp2 = tmp2 + 1;
            end
        end
        l1 = -l1/2;
        l2 = -l2/2;
        % log[i_1,i_2ごとの最尤推定量]
        L1 = Ri(1)*log(alpha) - Ri(1)*M*log(2*pi)/2 - Ri(1)*log(det(Vi1))/2 + l1;
        L2 = Ri(2)*log(alpha) - Ri(2)*M*log(2*pi)/2 - Ri(2)*log(det(Vi2))/2 + l2;
        % i_1,i_2を統一
        BIC = -2*(L1+L2) + M*(M+3)*log(R);
    end
end