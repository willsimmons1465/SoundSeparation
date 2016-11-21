THRESHOLD = 1000;
%THRESHOLD = 1;
WINDOWMILLIS = 1000;
%MINFREQ = 100; Not actually a constant
WF = 1;
WA = 1;
WH = 20;
MAXSINS = 100;
NUMSOURCES = 2;
TRAJECTORYDELTAFMAX = 0.1;

[samples, framerate] = audioread('C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Sound Samples\Square Intervals\3st.wav');
%[samples, framerate] = audioread('C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Sound Samples\Trumpet Violin\TA3_Vfs4.wav');

monoChannel = samples(:,1);

sampleLength = size(monoChannel);

windowSize = framerate*WINDOWMILLIS/1000;

spectrum = spectrogram(monoChannel, windowSize, floor(windowSize/2), windowSize, framerate);

%spectrogram(monoChannel, windowSize, windowSize/2, windowSize, framerate);

thresholdedSpectrum = abs(spectrum);
thresholdedSpectrum(thresholdedSpectrum < THRESHOLD) = 0;

%imagesc(thresholdedSpectrum)
%colormap(gray)

spectrumSize = size(spectrum);
%peaks = zeros(MAXSINS, spectrumSize(2));
%peakLocations = zeros(MAXSINS, spectrumSize(2));
%[p, pl] = findpeaks(thresholdedSpectrum(:,1))
%for i = 1:spectrumSize(2)
%    [peaks(:,i), peakLocations(:,i)] = findpeaks(thresholdedSpectrum(:,i));
%end

%Extract sinusoids
[peaks, peakLocations] = findpeaks(thresholdedSpectrum(:,1));
sinAmps = peaks(:);
sinFreqs = peakLocations(:);
for i = 2:spectrumSize(2)
    [peaks, peakLocations] = findpeaks(thresholdedSpectrum(:,i));
    numPeaks = size(peaks(:));
    numSins = size(sinFreqs(:,i-1));
    if numPeaks(1) == 0
        sinAmps(:, i) = 0;
        sinFreqs(:, i) = 0;
    end
    
    for j = 1:numPeaks(1)
        peakLoc = peakLocations(j);
        peak = peaks(j);
        matches = find(abs(log(sinFreqs(:,i-1)/peakLoc)) < log(1+TRAJECTORYDELTAFMAX));
        if isempty(matches)
            sinAmps(numSins+1,i) = peak;
            sinFreqs(numSins+1,i) = peakLoc;
        else
            sinAmps(matches, i) = peak;
            sinFreqs(matches, i) = peakLoc;
        end
    end
end

%sinAmps = peaks;
%sinFreqs = peakLocations;

sinAmpsSize = size(sinAmps);
numSins = sinAmpsSize(1);

meanAmps = mean(sinAmps, 2);
meanFreqs = mean(sinFreqs, 2);
for i = 1:spectrumSize(2)
    normAmps(:,i) = sinAmps(:,i) ./ meanAmps;
    normFreqs(:,i) = sinFreqs(:,i) ./ meanFreqs;
end

MinFreq = min(meanFreqs);

df = zeros(numSins, numSins);
da = zeros(numSins, numSins);
dh = Inf(numSins, numSins);

for i = 1:numSins
    for j = 1:numSins
        df(i,j) = sum((normFreqs(i,:) - normFreqs(j,:)).^2);
        da(i,j) = sum((normAmps(i,:) - normAmps(j,:)).^2);
        
        for a = 1:floor(meanFreqs(i)*(1+TRAJECTORYDELTAFMAX)/MinFreq)
            for b = 1:floor(meanFreqs(j)*(1+TRAJECTORYDELTAFMAX)/MinFreq)
                dhab = abs(log((meanFreqs(i)*b)/(meanFreqs(j)*a)));
                if dhab < dh(i,j)
                    dh(i,j) = dhab;
                end
            end
        end
    end
end

d = WF*df + WA*da + WH*dh;

%Identify Clusters
% minSplit = zeros(numSins,1);
% minD = Inf(1);
% for i = 1:NUMSOURCES^numSins
%     for j = 0:numSins-1
%         membershipVector(j+1) = mod(floor(i/NUMSOURCES^j), NUMSOURCES);
%     end
%     %membershipVector = de2bi(i, [], NUMSOURCES);
%     subsetDs = zeros(numSins, numSins);
%     
%     for j = 0:NUMSOURCES-1
%         %subsetDsIncr = d;
%         %subsetDsIncr(membershipVector ~= j, :) = 0;
%         mask = find(membershipVector == j);
%         subsetDs(mask, mask) = d(mask, mask);
%     end
%     
%     dSum = sum(subsetDs(:));
%     if (dSum < minD)
%         minSplit = membershipVector;
%         minD = dSum;
%     end
% end
minSplit = HierClusterN(d, NUMSOURCES);

%Synthesis
amps = zeros(numSins, sampleLength(1));
freqs = zeros(numSins, sampleLength(1));
expSinusoids = zeros(numSins, sampleLength(1));
outputSinusoids = zeros(numSins, sampleLength(1));

phases = ones(numSins, 1);
for i = 1:numSins
    nonZeroFreqs = sinFreqs(i, find(sinFreqs(i,:) ~= 0));
    phases(i) = spectrum(nonZeroFreqs(1));
end

%Interpolate amplitudes and frequencies
for i = 1:sampleLength(1)
    posInFrame = mod(i, windowSize);
    interpolationProportion = posInFrame / windowSize;
    amps(:, i) = interpolationProportion * sinAmps(:, floor(i/windowSize)+1) + (1-interpolationProportion) * sinAmps(:, ceil(i/windowSize)+1);
    freqs(:, i) = interpolationProportion * sinFreqs(:, floor(i/windowSize)+1) + (1-interpolationProportion) * sinFreqs(:, ceil(i/windowSize)+1);
    expSinusoids(:, i) = phases .* exp(2* pi * 1j * freqs(:,i) * i / framerate);
    outputSinusoids(:, i) = amps(:,i) .* real(expSinusoids(:,i));
end

outputSamples = zeros(NUMSOURCES, sampleLength(1));
for i = 0:NUMSOURCES-1
    mask = find(minSplit == i);
    outputSamples(i+1, :) = sum(outputSinusoids(mask, :));
end

peakAmps = max(outputSamples);
scaledOutputSamples = outputSamples / max(peakAmps);
audiowrite('C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Output Sounds SMProto\Square Intervals\3st0.wav', scaledOutputSamples(1, :), framerate);
audiowrite('C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Output Sounds SMProto\Square Intervals\3st1.wav', scaledOutputSamples(2, :), framerate);

%audiowrite('C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Output Sounds SMProto\Trumpet Violin\TA3_VFs4_0.wav', scaledOutputSamples(1, :), framerate);
%audiowrite('C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Output Sounds SMProto\Trumpet Violin\TA3_VFs4_1.wav', scaledOutputSamples(2, :), framerate);