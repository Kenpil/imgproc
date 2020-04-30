clear all;
clc;


width = 2737;
height = 1827;
load('Nachsholim/image_set_11/distanceFromCamera.mat');
depthMap = dist_map_r;
TF = isnan(dist_map_r);
depthMap(TF) = 100;
save('Nachsholim/image_set_11/depthMap5.mat', 'depthMap');

load('Nachsholim/image_set_11/depthMap5.mat');

figure;
% h = surf(-depthMap);
% zlim([-10 0])
depthlim = [0 15];
imagesc(depthMap,depthlim);
%set(h,'linestyle','none');