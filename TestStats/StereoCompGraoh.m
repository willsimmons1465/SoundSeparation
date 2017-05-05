stereo = [0, 0.2, 0.4, 0.6, 0.8, 1];
stereo = 360*asin(sqrt((1+stereo)/2))/pi - 90;
%sinSDR = [-19.3531, -19.3625, -19.4449, -3.45547, -2.68074, -3.29317];
%sinSIR = [5.21472, 5.28647, 5.49725, 27.5853, 24.0712, 25.7475];
%sinSAR = [-9.36108, -9.35276, -9.3573, -3.16943, -2.34064, -3.00836];
SDR = [-5.45978, -5.07527, -4.98738, -4.70931, -4.90512, -6.11516];
SDREr = [10.2959, 10.3967, 10.6855, 10.9995, 11.8473, 16.0631];
SIR = [31.4143, 31.3354, 32.219, 36.1511, 37.6237, 115.114];
SIREr = [22.1892, 22.2345, 22.8461, 24.0492, 25.8629, 255.13];
SAR = [-0.719761, -0.403282, -0.287791, -0.0166364, 0.228215, 0.29122];
SAREr = [5.41102, 5.37988, 5.56264, 5.64512, 5.52657, 5.39498];


figure
subplot(3, 1, 1)
%plot(stereo, sinSDR)
%errorbar(stereo, SDR, SDREr)
plot(stereo, SDR)
ylabel('SDR/dB')
%title('Separation performance against stereo separation of individual sounds')

subplot(3, 1, 2)
%plot(stereo, sinSIR)
%errorbar(stereo, SIR, SIREr)
plot(stereo, SIR)
%ylim([30, 45])
%xlim([0, 1])
ylim([30, 45])
xlim([0 90])
ylabel('SIR/dB')

subplot(3, 1, 3)
%plot(stereo, sinSAR)
%errorbar(stereo, SAR, SAREr)
plot(stereo, SAR)
ylabel('SAR/dB')
xlabel('Stereo position from centre/degrees')