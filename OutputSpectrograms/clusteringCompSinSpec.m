WINDOWSIZE = 1000;
HOPSIZE = 500;
NFFT = 2000;
FILE1 = 'C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Sound Samples\StandardTest\clarinet_As5_05_pianissimo_normal.wav';
FILE2 = 'C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Sound Samples\StandardTest\violin_C4_025_forte_arco-normal.wav';
FILE3 = 'C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\OutputSpectrograms\testOutputH_1171.wav';
FILE4 = 'C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\OutputSpectrograms\testOutputH_1170.wav';
FILE5 = 'C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\OutputSpectrograms\testOutput_1171.wav';
FILE6 = 'C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\OutputSpectrograms\testOutput_1170.wav';
FILE7 = 'C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\OutputSpectrograms\testOutputSM_1171.wav';
FILE8 = 'C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\OutputSpectrograms\testOutputSM_1170.wav';
FILE9 = 'C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\OutputSpectrograms\testOutputSN_1171.wav';
FILE10 = 'C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\OutputSpectrograms\testOutputSN_1170.wav';

[samples1, framerate] = audioread(FILE1);
[samples2, framerate] = audioread(FILE2);
[samples3, framerate] = audioread(FILE3);
[samples4, framerate] = audioread(FILE4);
[samples5, framerate] = audioread(FILE5);
[samples6, framerate] = audioread(FILE6);
[samples7, framerate] = audioread(FILE7);
[samples8, framerate] = audioread(FILE8);
[samples9, framerate] = audioread(FILE9);
[samples10, framerate] = audioread(FILE10);

if size(samples1, 2) == 2
    samples1 = samples1(:,1) + samples1(:,2);
end
if size(samples2, 2) == 2
    samples2 = samples2(:,1) + samples2(:,2);
end
if size(samples3, 2) == 2
    samples3 = samples3(:,1) + samples3(:,2);
end
if size(samples4, 2) == 2
    samples4 = samples4(:,1) + samples4(:,2);
end
if size(samples5, 2) == 2
    samples5 = samples5(:,1) + samples5(:,2);
end
if size(samples6, 2) == 2
    samples6 = samples6(:,1) + samples6(:,2);
end
if size(samples7, 2) == 2
    samples7 = samples7(:,1) + samples7(:,2);
end
if size(samples8, 2) == 2
    samples8 = samples8(:,1) + samples8(:,2);
end
if size(samples9, 2) == 2
    samples9 = samples9(:,1) + samples9(:,2);
end
if size(samples10, 2) == 2
    samples10 = samples10(:,1) + samples10(:,2);
end

length = max([size(samples1, 1), size(samples2, 1), size(samples3, 1), size(samples4, 1)]);
temp = samples1;
samples1 = zeros(length, 1);
samples1(1:size(temp, 1), :) = temp;
temp = samples2;
samples2 = zeros(length, 1);
samples2(1:size(temp, 1), :) = temp;
temp = samples3;
samples3 = zeros(length, 1);
samples3(1:size(temp, 1), :) = temp;
temp = samples4;
samples4 = zeros(length, 1);
samples4(1:size(temp, 1), :) = temp;
temp = samples5;
samples5 = zeros(length, 1);
samples5(1:size(temp, 1), :) = temp;
temp = samples6;
samples6 = zeros(length, 1);
samples6(1:size(temp, 1), :) = temp;
temp = samples7;
samples7 = zeros(length, 1);
samples7(1:size(temp, 1), :) = temp;
temp = samples8;
samples8 = zeros(length, 1);
samples8(1:size(temp, 1), :) = temp;
temp = samples9;
samples9 = zeros(length, 1);
samples9(1:size(temp, 1), :) = temp;
temp = samples10;
samples10 = zeros(length, 1);
samples10(1:size(temp, 1), :) = temp;

figure('Position', [100, -200, 600, 900])
subplot(3, 2, 1)
spectrogram(samples1, WINDOWSIZE, HOPSIZE, NFFT, framerate, 'yaxis', 'MinThreshold', -90);
cscale = caxis;
ylim([0 5]);

subplot(3, 2, 2)
spectrogram(samples2, WINDOWSIZE, HOPSIZE, NFFT, framerate, 'yaxis');
caxis(cscale);
ylim([0 5]);

subplot(3, 2, 3)
spectrogram(samples3, WINDOWSIZE, HOPSIZE, NFFT, framerate, 'yaxis');
caxis(cscale);
ylim([0 5]);

subplot(3, 2, 4)
spectrogram(samples4, WINDOWSIZE, HOPSIZE, NFFT, framerate, 'yaxis');
caxis(cscale);
ylim([0 5]);

subplot(3, 2, 5)
spectrogram(samples5, WINDOWSIZE, HOPSIZE, NFFT, framerate, 'yaxis');
caxis(cscale);
ylim([0 5]);

subplot(3, 2, 6)
spectrogram(samples6, WINDOWSIZE, HOPSIZE, NFFT, framerate, 'yaxis');
caxis(cscale);
ylim([0 5]);

cmap = [1, 1, 1
    1, 0.5, 0.5
    1, 0, 0
    0.5, 0, 0
    0, 0, 0];
colormap(cmap);

ch=findall(gcf,'tag','Colorbar');
delete(ch);

figure('Position', [100, -200, 600, 435]);
subplot(2, 2, 1)
spectrogram(samples7, WINDOWSIZE, HOPSIZE, NFFT, framerate, 'yaxis');
caxis(cscale);
ylim([0 5]);

subplot(2, 2, 2)
spectrogram(samples8, WINDOWSIZE, HOPSIZE, NFFT, framerate, 'yaxis');
caxis(cscale);
ylim([0 5]);

subplot(2, 2, 3)
spectrogram(samples9, WINDOWSIZE, HOPSIZE, NFFT, framerate, 'yaxis');
caxis(cscale);
ylim([0 5]);

subplot(2, 2, 4)
spectrogram(samples10, WINDOWSIZE, HOPSIZE, NFFT, framerate, 'yaxis');
caxis(cscale);
ylim([0 5]);

cmap = [1, 1, 1
    1, 0.5, 0.5
    1, 0, 0
    0.5, 0, 0
    0, 0, 0];
colormap(cmap);

ch=findall(gcf,'tag','Colorbar');
delete(ch);