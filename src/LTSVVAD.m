function [ postVAD, preVAD ] = LTSVVAD(s, fs, wsec, isec, parameter, enhance, thr)
% set the default threshold if none supplied
if(nargin < 7)
    thr = 10;
end

% default - no speech enhancement
if(nargin < 6)
    enhance = 0;
end

if(nargin < 5)
    M = 2;
    R = 6;
else
    M = parameter.M;
    R = parameter.R;
end
% -----------------------------------------------------------
% PARAMETERS
% -----------------------------------------------------------
nDFT = 2048; % DFT number
K1 = round(nDFT * 500 / fs);
K2 = round(nDFT * 4000 / fs);
if isec == 0.01
    B = 14; % buffer length
    Sp = 3; % speech possible
    Sl = 5; % speech likely
    Ls = 18; % short hangover time
    Lm = 25; % medium hangover time
else
    B = 7; % buffer length
    Sp = 2; % speech possible
    Sl = 3; % speech likely
    Ls = 5; % short hangover time
    Lm = 8; % medium hangover time
end

% -----------------------------------------------------------
% PRE-PROCESSING
% -----------------------------------------------------------
% perform speech enhancement if necessary
if(enhance == 1)
    s = specsub(s,fs);
end
if(enhance == 2)
    s = ssubmmse(s,fs);
end

% number of samples per window
winSamples = round(wsec*fs);
incSamples = round(isec*fs);
% enframe the signal using hamming window
frames = enframe(s,hamming(winSamples,'periodic'),incSamples);

% -----------------------------------------------------------
% FEATURE EXTRACTION
% -----------------------------------------------------------
% calculate the spectrum for each frame
dft = rfft(frames,nDFT,2);
% calculate the amplitude spectrum for each frame
amplitudeSpectrum = abs(dft);
% calculate the Power Spectrum of the noisy signal
signalPS = amplitudeSpectrum .^ 2;
% estimate the signal spectrum by averaging spectral estimates 
% of M consecutive frames
noFrames = size(signalPS,1);
specBW = zeros(noFrames, K2 - K1 + 1);
for n = 1 : M - 1
    specBW(n, :) = sum(signalPS(1 : n, K1 : K2), 1) / n;
end
for n = M : noFrames
    specBW(n, :) = sum(signalPS((n - M + 1) : n, K1 : K2), 1) / M;
end
% entropy measure on the normalized short-time spectrum computed
 %at frequency w_k over R consecutive frames
entropy = zeros(noFrames, K2 - K1 + 1);
for m = 1 : R - 1
    sumSpec = repmat(sum(specBW(1 : m, :), 1), m, 1);
    ratioSpec = specBW(1 : m, :) ./ sumSpec;
    entropy(m, :) = - sum(ratioSpec .* log(ratioSpec));
end
for m = R : noFrames
    sumSpec = repmat(sum(specBW((m - R + 1) : m, :), 1), R, 1);
    ratioSpec = specBW((m - R + 1) : m, :) ./ sumSpec;
    entropy(m, :) = - sum(ratioSpec .* log(ratioSpec));
end
% The long-term signal variability (LTSV) measure
LTSV = zeros(noFrames, 1);
for m = 1 : noFrames
    LTSV(m) = log10(sum((entropy(m, :) - mean(entropy(m, :))).^2) / (K2 - K1 + 1));
end
% -----------------------------------------------------------
% CLASSIFICATION
% -----------------------------------------------------------
thrno = length(thr);
framesVAD = zeros(thrno,length(LTSV));
preVAD = zeros(thrno,length(s));
for j = 1 : thrno
    % calculate the VAD decisions
    framesVAD(j,LTSV > thr(j)) = 1;
    framesVAD(j,LTSV <= thr(j)) = 0;
    preVAD(j,1:noFrames*incSamples) = reshape(repmat(framesVAD(j,:),...
        incSamples, 1), 1, noFrames*incSamples);
end

% -----------------------------------------------------------
% POST-PROCESSING
% -----------------------------------------------------------
% apply hang-over scheme from the original paper
hangoverVAD = zeros(thrno,noFrames);
postVAD = zeros(thrno,length(s));
for j = 1 : thrno
    T = 0;
    for i = 1:(noFrames-B+1)
        M = maxConsOnes(framesVAD(j,i:i+B-1));

        if(M >= Sl)
            T = Lm;
        elseif(M >= Sp && T < Ls)
            T = Ls;
        elseif(M < Sp && T > 0)
            T = T - 1;
        end

        if(T > 0)
            hangoverVAD(j,i) = 1;
        end
    end
    hangoverVAD(j,noFrames-B+2:noFrames) = framesVAD(j,noFrames-B+2:noFrames);
    % transform the VAD frames to samples
    postVAD(j,1:noFrames*incSamples) = reshape(repmat(hangoverVAD(j,:),...
        incSamples, 1), 1, noFrames*incSamples);
end
end

function max = maxConsOnes(seq)
    M = 0;
    max = 0;
    
    for i = 1:length(seq)
        if(seq(i) == 1)
            M = M+1;
            if(M > max)
                max = M;
            end
        else
            M = 0;
        end
    end
end