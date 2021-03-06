categories = {'Without (Sin)', 'With (Sin)', 'Without (NMF)', 'With (NMF)'};
SDRs = [-34.778, -11.3366, -8.94168, -5.90039;
    -9.08361, -0.0256131, -1.97583, 0.338262;
    -2.09438, 1.81522, 0.30395, 2.61904;
    1.7839, 4.04335, 5.70751, 5.83092;
    10.5181, 13.2481, 12.4805, 13.7146];
SIRs = [-8.49768, 0.0151946, -1.1504, -1.10483;
    14.9408, 4.79472, 1.72184, 2.54299;
    30.2107, 12.1206, 7.06244, 8.07935;
    45.3237, 19.3327, 18.822, 12.5642;
    93.9436, 48.5012, 34.1076, 24.0849];
SARs = [-10.3683, 0.540008, 0.104375, 2.04667;
    -3.64178, 3.25905, 3.91289, 4.87445;
    -0.452829, 4.91571, 5.20714, 8.60272;
    3.66214, 6.62454, 9.58017, 13.9642;
    12.281, 35.6743, 12.602, 23.4121];

figure
subplot(3, 1, 1)
boxplot(SDRs, categories)
ylabel('SDR/dB')
xticks([]);
%title('Separation performance with and without the reversibility property')

subplot(3, 1, 2)
boxplot(SIRs, categories)
ylabel('SIR/dB')
xticks([]);

subplot(3, 1, 3)
boxplot(SARs, categories, 'Whisker', 2)
ylabel('SAR/dB')