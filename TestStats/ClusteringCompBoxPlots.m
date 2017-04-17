cOp = {'Hard', 'Soft', 'NMF', 'Naive'};
sinSDR = [-27.0381, -29.3136, -26.1243, -27.8911;
    -6.47037, -6.95663, -7.88166, -8.27079;
    -1.77345, -1.23999, -3.84938, -3.00128;
    1.7084, 1.61314, 0.467746, 1.02411;
    10.508, 10.508, 7.12542, 8.7251];
sinSIR = [3.35525, 4.64453, 0.379755, 0.508087;
    23.1529, 23.4161, 4.02963, 5.2152;
    31.2021, 31.1084, 6.36298, 12.4585;
    54.1427, 57.1629, 15.9052, 26.1536;
    91.9698, 91.9698, 50.6097, 70.4597];
sinSAR = [-9.95423, -9.85423, -4.74493, -7.68576;
    -3.42425, -2.29602, -0.778708, -2.25727;
    -1.19609, -0.428729, 2.58864, -1.4044;
    1.93284, 2.34062, 4.30558, 1.8299;
    10.5081, 10.5081, 8.92478, 8.72512];
matSDR = [-11.7677, -8.94168, -14.3621;
    -3.48341, -1.97583, -4.24805;
    2.18627, 0.30395, -0.305671;
    7.31409, 5.70751, 3.73385;
    12.4822, 12.4805, 11.2295];
matSIR = [1.90226, -1.1504, -1.35573;
    7.38514, 1.72184, 4.7028;
    20.8073, 7.06244, 10.8948;
    26.2854, 18.822, 17.868;
    52.5768, 34.1076, 51.8927];
matSAR = [-7.67816, 0.104375, -12.9222;
    1.26376, 3.91289, 1.15286;
    4.1514, 5.20714, 3.98315;
    8.43368, 9.58017, 8.58706;
    12.506, 12.602, 15.5723];

cOpAll = {'Hard (S)', 'Soft (S)', 'Matrix (S)', 'Naive (S)', 'Hard (MF)', 'Soft (MF)', 'Matrix (MF)'};
SDR = [sinSDR, matSDR];
SIR = [sinSIR, matSIR];
SAR = [sinSAR, matSAR];

figure
subplot(3, 1, 1)
%boxplot(sinSDR, cOp)
%hold on
%boxplot(matSDR, cOp(1:3))
boxplot(SDR, cOpAll)
ylabel('SDR')
xticks([])
%legend('Sinusoids', 'Matrix Factors')
title('Separation performance against clustering method')

subplot(3, 1, 2)
%boxplot(sinSIR, cOp)
%hold on
%boxplot(matSIR, cOp(1:3))
boxplot(SIR, cOpAll)
ylabel('SIR')
xticks([])

subplot(3, 1, 3)
%boxplot(sinSAR, cOp)
%hold on
%boxplot(matSAR, cOp(1:3))
boxplot(SAR, cOpAll)
xlabel('Clustering Method')
ylabel('SAR')