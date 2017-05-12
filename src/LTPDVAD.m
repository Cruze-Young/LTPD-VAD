function [ postVAD, preVAD ] = LTPDVAD(s, fs, wsec, isec, parameter, enhance, thr)
% set the default threshold if none supplied
if(nargin < 7)
    thr = 10;
end

% default - no speech enhancement
if(nargin < 6)
    enhance = 0;
end

if(nargin < 5)
    N = 3;
else
    N = parameter.N;
end

% -----------------------------------------------------------
% PARAMETERS
% -----------------------------------------------------------
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

incSamples = round(isec*fs);
% -----------------------------------------------------------
% FEATURE EXTRACTION
% -----------------------------------------------------------
paramAudio.useResampling = 1;
[f_audio,sideinfo] = waveform_to_audio(s,fs,16, paramAudio);
shiftFB = estimateTuning(f_audio);

paramPitch.winLenSTMSP = round(wsec*22050);
paramPitch.winOvSTMSP = round((wsec-isec)*22050);
paramPitch.shiftFB = shiftFB;
[f_pitch,~] = audio_to_pitch_via_FB(f_audio,paramPitch,sideinfo);
    
if size(f_pitch,2)*incSamples > length(s)
    f_pitch = f_pitch(:, 1:end-1);
    if size(f_pitch,2)*incSamples > length(s)
        f_pitch = f_pitch(:, 2:end);
    end
end
f_pitch = f_pitch(67:102,:)';
% estimate the Power Spectrum of the noise
noisePS = estnoiseg(f_pitch,isec);
avgNoisePS = mean(noisePS)+eps;

noFrames = size(f_pitch,1);
freqBins = size(f_pitch,2);
LTSE = zeros(noFrames,freqBins);
LTPD = zeros(noFrames,1);

for i = 1:noFrames
    for j = 1:freqBins
        startIndex = max(i-N,1);
        endIndex = min(i+N,noFrames);
        
        LTSE(i,j) = max(f_pitch(startIndex:endIndex,j));
    end
    
    currentLTSEsq = LTSE(i,:).^2;
    LTSEnoise = currentLTSEsq./avgNoisePS;
    LTPD(i) = 10*log10(mean(LTSEnoise));
end

% -----------------------------------------------------------
% CLASSIFICATION
% -----------------------------------------------------------
thrno = length(thr);
framesVAD = zeros(thrno,length(LTPD));
preVAD = zeros(thrno,length(s));
for j = 1 : thrno
    % calculate the VAD decisions
    framesVAD(j,LTPD > thr(j)) = 1;
    framesVAD(j,LTPD <= thr(j)) = 0;
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