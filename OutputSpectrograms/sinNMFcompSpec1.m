WINDOWSIZE = 1000;
HOPSIZE = 500;
NFFT = 2000;
FILE1 = 'C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Sound Samples\StandardTest\clarinet_B6_1_pianissimo_normal.wav';
FILE2 = 'C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Sound Samples\StandardTest\violin_Ds5_1_mezzo-forte_arco-normal.wav';
FILE3 = 'C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\OutputSpectrograms\testOutput_2181.wav';
FILE4 = 'C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\OutputSpectrograms\testOutput_2180.wav';
FILE5 = 'C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\OutputSpectrograms\testOutputM_2181.wav';
FILE6 = 'C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\OutputSpectrograms\testOutputM_2180.wav';

[samples1, framerate] = audioread(FILE1);
[samples2, framerate] = audioread(FILE2);
[samples3, framerate] = audioread(FILE3);
[samples4, framerate] = audioread(FILE4);
[samples5, framerate] = audioread(FILE5);
[samples6, framerate] = audioread(FILE6);

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