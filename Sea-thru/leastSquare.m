clear all;
clc;

x = (0:0.01:10)';
BdataFull = 0.4*(1-exp(-2.2*x))+0.55*exp(-3.6*x);
bEstimate = fittype(@(A,B,C,D,x) A*(1-exp(-1*B*x))+C*exp(-1*D*x));
[Bfit, gof] = fit(x, BdataFull, bEstimate, 'StartPoint', [0.5, 2.5, 0.5, 3.0]);
figure;
plot(x,BdataFull);
hold on;
plot(x,Bfit(x));
hold off;