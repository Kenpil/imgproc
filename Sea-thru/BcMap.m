clear all;
clc;


width = 2737;
height = 1827;
load('Nachsholim/image_set_01/distanceFromCamera.mat');
depthMap = dist_map_r;
TF = isnan(dist_map_r);
depthMap(TF) = 100;
save('Nachsholim/image_set_01/depthMap.mat', 'depthMap');

load('Nachsholim/image_set_01/depthMap.mat');

figure;
h = surf(-depthMap);
zlim([-10 0])
set(h,'linestyle','none');