%THRESHOLD = 200;
THRESHOLD = 30;
WINDOWMILLIS = 32768*1000/44100;
%WINDOWMILLIS = 1000;
%MINFREQ = 100; Not actually a constant
WF = 1;
WA = 2.5;
WH = 1000;
%MAXSINS = 100; Unused
NUMSOURCES = 2;
TRAJECTORYDELTAFMAX = 0.02;
WINDOWOVERLAP = 0.5;

%[samples, framerate] = audioread('C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Sound Samples\Square Intervals\3st.wav');
[samples, framerate] = audioread('C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Sound Samples\Trumpet Violin\TA3_Vfs4.wav');

monoChannel = samples(:,1);

sampleLength = size(monoChannel);

windowSize = framerate*WINDOWMILLIS/1000;
overlap = floor(windowSize*WINDOWOVERLAP);

%spectrum = spectrogram(monoChannel, windowSize, floor(windowSize/2), windowSize, framerate);
spectrum = spectrogram(monoChannel, windowSize, overlap, framerate, framerate);


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
sinFreqBins = peakLocations(:);
for i = 2:spectrumSize(2)
    [peaks, peakLocations] = findpeaks(thresholdedSpectrum(:,i));
    numPeaks = size(peaks(:));
    numSins = size(sinFreqBins(:,i-1));
    if numPeaks(1) == 0
        sinAmps(:, i) = 0;
        sinFreqBins(:, i) = 0;
    end
    
    for j = 1:numPeaks(1)
        peakLoc = peakLocations(j);
        peak = peaks(j);
        matches = find(abs(log(sinFreqBins(:,i-1)/peakLoc)) < log(1+TRAJECTORYDELTAFMAX));
        if isempty(matches)
            sinAmps(numSins(1)+1,i) = peak;
            sinFreqBins(numSins(1)+1,i) = peakLoc;
            numSins = size(sinFreqBins(:,i-1));
        else
            sinAmps(matches, i) = peak;
            sinFreqBins(matches, i) = peakLoc;
        end
    end
end
sinFreqs = sinFreqBins * framerate / framerate;

%sinAmps = peaks;
%sinFreqs = peakLocations;

sinAmpsSize = size(sinAmps);
numSins = sinAmpsSize(1);

meanAmps = MeanIgnoringZeros(sinAmps);
meanFreqs = MeanIgnoringZeros(sinFreqs);
normAmps = zeros(sinAmpsSize);
normFreqs = zeros(sinAmpsSize);
for i = 1:spectrumSize(2)
    normAmps(:,i) = sinAmps(:,i) ./ meanAmps;
    normFreqs(:,i) = sinFreqs(:,i) ./ meanFreqs;
end

MinFreq = min(meanFreqs);

df = zeros(numSins, numSins);
da = zeros(numSins, numSins);
dh = Inf(numSins, numSins);

for i = 1:numSins
    for j = i:numSins
        df(i,j) = sum((normFreqs(i,:) - normFreqs(j,:)).^2);
        da(i,j) = sum((normAmps(i,:) - normAmps(j,:)).^2);
        
        maxB = floor(meanFreqs(j)*(1+TRAJECTORYDELTAFMAX)/MinFreq);
        
        for a = 1:floor(meanFreqs(i)*(1+TRAJECTORYDELTAFMAX)/MinFreq)
%             for b = 1:floor(meanFreqs(j)*(1+TRAJECTORYDELTAFMAX)/MinFreq)
%                 dhab = abs(log((meanFreqs(i)*b)/(meanFreqs(j)*a)));
%                 if dhab < dh(i,j)
%                     dh(i,j) = dhab;
%                 end
%             end
            b = a*meanFreqs(j)/meanFreqs(i);
            if b >= maxB+1
                break;
            end
            
            if b >= 1
                dhab = abs(log((meanFreqs(i)*floor(b))/(meanFreqs(j)*a)));
                if dhab < dh(i,j)
                    dh(i,j) = dhab;
                end
            end
            
            if b > maxB
                break;
            end
            
            dhab = abs(log((meanFreqs(i)*ceil(b))/(meanFreqs(j)*a)));
            if dhab < dh(i,j)
                dh(i,j) = dhab;
            end
        end
        
        df(j,i) = df(i,j);
        da(j,i) = da(i,j);
        dh(j,i) = dh(i,j);
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
    nonZeroFreqs = sinFreqBins(i, sinFreqBins(i,:) ~= 0);
    phases(i) = spectrum(nonZeroFreqs(1));
    phases(i) = phases(i)/abs(phases(i));
end

%Interpolate amplitudes and frequencies
% for i = 1:sampleLength(1)
%     posInFrame = mod(i, windowSize);
%     interpolationProportion = posInFrame / windowSize;
%     amps(:, i) = interpolationProportion * sinAmps(:, floor(i/windowSize)+1) + (1-interpolationProportion) * sinAmps(:, ceil(i/windowSize)+1);
%     freqs(:, i) = interpolationProportion * sinFreqs(:, floor(i/windowSize)+1) + (1-interpolationProportion) * sinFreqs(:, ceil(i/windowSize)+1);
%     expSinusoids(:, i) = phases .* exp(2* pi * 1j * freqs(:,i) * i / framerate);
%     outputSinusoids(:, i) = amps(:,i) .* real(expSinusoids(:,i));
% end
sinAmpsExtended = [zeros(numSins,1) sinAmps zeros(numSins,1)];
sinFreqsExtended = [sinFreqs(:,1) sinFreqs sinFreqs(:,end)];
for i = 1:sampleLength(1)
    posInFrame = mod(i, windowSize-overlap);
    baseFrame = floor(i/(windowSize-overlap))+1;
    upperFrame = ceil(i/(windowSize-overlap))+1;
    if upperFrame > size(sinAmpsExtended, 2)
        baseFrame = size(sinAmpsExtended, 2);
        upperFrame = size(sinAmpsExtended, 2);
    end
    interpolationProportion = posInFrame / (windowSize-overlap);
    amps = (1-interpolationProportion) * sinAmpsExtended(:, baseFrame) + interpolationProportion * sinAmpsExtended(:, upperFrame);
    freqs = (1-interpolationProportion) * sinFreqsExtended(:, baseFrame) + interpolationProportion * sinFreqsExtended(:, upperFrame);
    expSinusoids = phases .* exp(2* pi * 1j * freqs * i / framerate);
    outputSinusoids(:, i) = amps .* real(expSinusoids);
end

outputSamples = zeros(NUMSOURCES, sampleLength(1));
for i = 0:NUMSOURCES-1
    mask = find(minSplit == i);
    outputSamples(i+1, :) = sum(outputSinusoids(mask, :));
end

peakAmps = max(outputSamples);
scaledOutputSamples = outputSamples / max(peakAmps);
%audiowrite('C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Output Sounds SMProto\Square Intervals\3st0.wav', scaledOutputSamples(1, :), framerate);
%audiowrite('C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Output Sounds SMProto\Square Intervals\3st1.wav', scaledOutputSamples(2, :), framerate);

audiowrite('C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Output Sounds SMProto\Trumpet Violin\TA3_VFs4_0.wav', scaledOutputSamples(1, :), framerate);
audiowrite('C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Output Sounds SMProto\Trumpet Violin\TA3_VFs4_1.wav', scaledOutputSamples(2, :), framerate);