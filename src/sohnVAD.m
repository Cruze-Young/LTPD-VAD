function [ postVAD, preVAD ] = sohnVAD(s,fs,wsec,isec,parameter,enhance,thr)
% set the default threshold if none supplied
if(nargin < 6)
    thr = 2;
end

% default - no speech enhancement
if(nargin < 5)
    enhance = 0;
end

% -----------------------------------------------------------
% PARAMETERS
% -----------------------------------------------------------
alpha = 0.95;
B = 7; % buffer length
Sp = 2; % speech possible
Sl = 3; % speech likely
Ls = 5; % short hangover time
Lm = 8; % medium hangover time

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
% enframe the signal using hanning window
frames = enframe(s,hanning(winSamples,'periodic'),incSamples);
noFrames = size(frames,1);

% -----------------------------------------------------------
% FEATURE EXTRACTION
% -----------------------------------------------------------
% calculate the DFT of each frame
dft = rfft(frames,winSamples,2);
% calculate the Power Spectrum of the noisy signal
signalPS = dft.*conj(dft);
% estimate the Power Spectrum of the noise
noisePS = estnoiseg(signalPS,wsec);
% calculate the a posteriori SNR
gamma = signalPS./noisePS;

A = 1;
noi = 1;
logRatio = zeros(noFrames,1);
for i = 1:noFrames
    gami = gamma(i,:);
    xi = alpha*((A.^2)./noi) + (1-alpha).*max(gami-1,0);
    V = (gami.*xi)./(1+xi);
    V2 = V/2;
    A = 0.5*sqrt(pi)*sqrt(V)./gami.*exp(-V2).*((1+V).*besseli(0,V2) + V.*besseli(1,V2));
    A(isnan(A)|isinf(A)) = 1;
    A = A.*sqrt(signalPS(i,:));
    noi = noisePS(i,:);
    
    logRatio(i) = mean(V - log(1+xi));
end

% -----------------------------------------------------------
% CLASSIFICATION
% -----------------------------------------------------------
% set the threshold
%thr = 2;
thrno = length(thr);
framesVAD = zeros(thrno,length(logRatio));
preVAD = zeros(thrno,length(s));
for j = 1 : thrno
    % calculate the VAD decisions
    framesVAD(j,logRatio > thr(j)) = 1;
    framesVAD(j,logRatio <= thr(j)) = 0;
    preVAD(1:noFrames*incSamples) = reshape(repmat(framesVAD(j,:),...
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