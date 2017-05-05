offset = [0, 50, 100, 150, 200];
sinSDR = [-12.3959, -7.41124, -7.16866, -8.21264, -7.71807];
sinSIR = [17.6314, 16.3814, 17.3437, 19.0665, 20.2695];
sinSAR = [-11.7158, -4.53312, -4.62139, -5.57893, -5.11637];
matSDR = [-3.21057, -3.1904, -1.51928, 0.625772, -8.63883];
matSIR = [8.55148, 5.07654, 4.46536, 3.87947, 3.57002];
matSAR = [-0.438837, -0.182423, 1.72026, 5.53185, -1.17835];
%sinSDR = [-5.05571, -4.78194, -2.59203, -4.84077, -5.1609];
%sinSIR = [31.5278, 28.0949, 30.6242, 31.2684, 39.1434];
%sinSAR = [-3.86843, -3.58969, -1.65403, -3.68113, -3.99888];
%matSDR = [0.0330537, 0.664474, 1.37812, 0.703439, 0.830332];
%matSIR = [12.9222, 15.5611, 15.4831, 16.7703, 10.7918];
%matSAR = [1.36485, 1.75665, 2.31825, 1.36264, 2.63884];


figure
subplot(3, 1, 1)
plot(offset, sinSDR)
hold on
plot(offset, matSDR)
ylabel('SDR/dB')
xticks([])
ylim([-20 10])
legend('Sin', 'NMF')
%title('Separation performance against relative offset of sounds')

subplot(3, 1, 2)
plot(offset, sinSIR)
hold on
plot(offset, matSIR)
ylabel('SIR/dB')
xticks([])
ylim([0 30])

subplot(3, 1, 3)
plot(offset, sinSAR)
hold on
plot(offset, matSAR)
ylabel('SAR/dB')
xlabel('Relative offset of sounds/ms')
ylim([-20 10])