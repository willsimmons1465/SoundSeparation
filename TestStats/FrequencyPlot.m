freq = [0, 10, 20, 30, 40, 50];
sinSDR = [11.6676, -22.0517, -13.0079, 11.7505, 12.8779, 14.0428];
sinSIR = [13.1501, -4.27103, 7.80673, 31.6438, 39.8788, 51.8908];
sinSAR = [17.5439, -2.13721, -2.83289, 11.8272, 12.9008, 14.0622];

matSDR = [13.3236, 0.00725665, 0.00129434, -0.00424634, -0.00789903, -0.00599553];
matSIR = [13.3735, 0.0351281, 0.0087878, 9.2862e-5, 0.00187707, 0.00141032];
matSAR = [34.0998, 27.1559, 32.8815, 33.5165, 31.7217, 32.9378];

figure
subplot(3, 1, 1)
plot(freq, sinSDR)
hold on
plot(freq, matSDR)
ylabel('SDR')
xticks([])
ylim([-30 20])
legend('Sinusoids', 'Matrix Factors')
title('Separation performance against relative frequency')

subplot(3, 1, 2)
plot(freq, sinSIR)
hold on
plot(freq, matSIR)
ylabel('SIR')
xticks([])
ylim([-20 60])

subplot(3, 1, 3)
plot(freq, sinSAR)
hold on
plot(freq, matSAR)
ylabel('SAR')
xlabel('Difference between base frequencies/Hz')
ylim([-10 40])