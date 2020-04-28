clear all;
clc;


width = 2737;
height = 1827;
load('Nachsholim/image_set_05/distanceFromCamera.mat');
depthMap = dist_map_l;
TF = isnan(dist_map_l);
depthMap(TF) = 100;
save('Nachsholim/image_set_05/depthMap3.mat', 'depthMap');

load('Nachsholim/image_set_05/depthMap3.mat');

figure;
% h = surf(-depthMap);
% zlim([-10 0])
depthlim = [0 15];
imagesc(depthMap,depthlim);
%set(h,'linestyle','none');