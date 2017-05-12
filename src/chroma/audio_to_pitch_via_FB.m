function [f_pitch,sideinfo] = audio_to_pitch_via_FB(f_audio,parameter,sideinfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: audio_to_pitch_via_FB
% Date of Revision: 2011-03
% Programmer: Meinard Mueller, Sebastian Ewert
%
% Description:
% Computing and saving of pitch features via a pre-designed filterbank. 
% features. For each window length specified via parameter.winLenSTMSP 
% the following is computed:
% - STMSP (short-time mean-square power) for each MIDI pitch between
%   parameter.midiMin and parameter.midiMax
% - STMSP subbands are stored in f_pitch, where f_pitch(p,:) contains
%   STMSP of subband of pitch p
% - sideinfo contains information of original pcm, which is saved along
%   with f_pitch into a single mat-file
% - Information f_pitch and sideinfo is stored in mat-file:
%   save(strcat(parameter.saveDir,parameter.saveFilename),'f_pitch','sideinfo');
%
% Input:
%        f_audio
%        parameter.winLenSTMSP = 4410;
%        parameter.shiftFB = 0;
%        parameter.midiMin = 21;
%        parameter.midiMax = 108;
%        parameter.save = 0;
%           parameter.saveDir = '';
%           parameter.saveFilename = '';
%           parameter.saveAsTuned = 0;
%        parameter.fs = 22050;
%        parameter.visualize = 0;
%
% Required files:
%        'MIDI_FB_ellip_pitch_60_96_22050_Q25.mat'
%        'MIDI_FB_ellip_pitch_60_96_22050_Q25_minusHalf.mat'
%        'MIDI_FB_ellip_pitch_60_96_22050_Q25_minusQuarter.mat'
%        'MIDI_FB_ellip_pitch_60_96_22050_Q25_minusThird.mat'
%        'MIDI_FB_ellip_pitch_60_96_22050_Q25_minusThreeQuarters.mat'
%        'MIDI_FB_ellip_pitch_60_96_22050_Q25_minusTwoThird.mat'
%
% Output:
%        f_pitch
%        sideinfo
%
%
% License:
%     This file is part of 'Chroma Toolbox'.
% 
%     'Chroma Toolbox' is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 2 of the License, or
%     (at your option) any later version.
% 
%     'Chroma Toolbox' is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with 'Chroma Toolbox'. If not, see <http://www.gnu.org/licenses/>.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
centerFreq = [   8.66;    9.18;    9.72;   10.30;   10.91;   11.56;   12.25;   12.98;...
                13.75;   14.57;   15.43;   16.35;   17.32;   18.35;   19.45;   20.60;...
                21.83;   23.12;   24.50;   25.96;   27.50;   29.14;   30.87;   32.70;...
                34.65;   35.71;   38.89;   41.20;   43.65;   46.25;   49.00;   51.91;...
                55.00;   56.27;   61.74;   65.41;   69.30;   73.42;   77.78;   82.41;...
                87.31;   92.50;   98.00;  103.83;  110.00;  116.54;  123.47;  130.81;...
               138.69;  146.83;  155.56;  164.81;  174.61;  185.00;  196.00;  207.65;...
               220.00;  233.08;  246.94;  261.63;  277.18;  293.66;  311.13;  329.63;...
               349.23;  369.99;  388.00;  415.30;  440.00;  466.16;  493.88;  523.25;...
               554.37;  567.33;  622.25;  659.26;  698.46;  739.99;  783.99;  830.61;...
               880.00;  932.33;  987.77; 1046.50; 1108.73; 1174.66; 1244.51; 1318.51;...
              1396.91; 1479.96; 1567.96; 1661.22; 1760.00; 1864.55; 1975.53; 2083.00;...
              2217.46; 2319.32; 2489.02; 2537.02; 2793.83; 2969.96; 3135.96; 3322.44;...
              3520.00; 3729.31; 3951.07; 4186.01; 4434.92; 4698.64; 4978.03; 5274.04;...
              5587.65; 5919.91; 6271.93; 6644.88; 7040.00; 7458.62; 7902.13; 8372.02;...
              8869.84; 9397.27; 9956.06; 10548.08; 11175.30; 11839.82; 12543.85] ;
if nargin<3
    sideinfo=[];
end

if nargin<2
    parameter=[];
end
if isfield(parameter,'visualize')==0
    parameter.visualize = 0;
end
if isfield(parameter,'save')==0
    parameter.save = 0;
end
if isfield(parameter,'saveDir')==0
    parameter.saveDir = '';
end
if isfield(parameter,'saveFilename')==0
    parameter.saveFilename = '';
end
if isfield(parameter,'saveAsTuned')==0
    parameter.saveAsTuned = 0;
end
if isfield(parameter,'fs')==0
    parameter.fs = 22050;
else
    if parameter.fs ~= 22050
        error('audio_to_pitch_via_FB not implemented yet for sample rates other than 22050.');
    end
end
if isfield(parameter,'midiMin')==0
    parameter.midiMin = 21;
end
if isfield(parameter,'midiMax')==0
    parameter.midiMax = 108;
%     parameter.midiMax = find(centerFreq > sideinfo.wav.sourceFs/2, 1, 'first');
end
if isfield(parameter,'winLenSTMSP')==0
    parameter.winLenSTMSP = 4410;
    %parameter.winLenSTMSP = [882 4410];
end
if isfield(parameter,'winOvSTMSP')==0
    parameter.winOvSTMSP = round(parameter.winLenSTMSP/2);
end
if isfield(parameter,'shiftFB')==0
    parameter.shiftFB = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parameter.shiftFB == 0
    load MIDI_FB_ellip_pitch_60_96_22050_Q25.mat
elseif parameter.shiftFB == 1
    load MIDI_FB_ellip_pitch_60_96_22050_Q25_minusQuarter.mat
elseif parameter.shiftFB == 2
    load MIDI_FB_ellip_pitch_60_96_22050_Q25_minusThird.mat
elseif parameter.shiftFB == 3
    load MIDI_FB_ellip_pitch_60_96_22050_Q25_minusHalf.mat
elseif parameter.shiftFB == 4
    load MIDI_FB_ellip_pitch_60_96_22050_Q25_minusTwoThird.mat
elseif parameter.shiftFB == 5
    load MIDI_FB_ellip_pitch_60_96_22050_Q25_minusThreeQuarters.mat
else
    error('Wrong shift parameter!')
end

fs_pitch = zeros(1,128);
fs_index = zeros(1,128);

fs_pitch(21:59) = 882;
fs_pitch(60:95) = 4410;
fs_pitch(96:120) = 22050;

fs_index(21:59) = 3;
fs_index(60:95) = 2;
fs_index(96:120) = 1;

pcm_ds = cell(3,1);
pcm_ds{1} = f_audio;
pcm_ds{2} = resample(pcm_ds{1},1,5,100);
pcm_ds{3} = resample(pcm_ds{2},1,5,100);

% fprintf('Computing subbands and STMSP for all pitches: (%i-%i): %4i',parameter.midiMin,parameter.midiMax,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute features for all pitches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

winLenSTMSP =  parameter.winLenSTMSP;
% winOvSTMSP  =  round(winLenSTMSP/2);
% winOvSTMSP = 0;
winOvSTMSP = parameter.winOvSTMSP;
featureRate =  parameter.fs./(winLenSTMSP-winOvSTMSP);  %formerly win_res
wav_size = size(f_audio,1);

num_window = length(winLenSTMSP);
f_pitch_energy = cell(num_window,1);
seg_pcm_num = cell(num_window,1);
seg_pcm_start = cell(num_window,1);
seg_pcm_stop = cell(num_window,1);
for w=1:num_window;
    step_size = winLenSTMSP(w)-winOvSTMSP(w);
    group_delay = round(winLenSTMSP(w)/2);
    seg_pcm_start{w} = [1 1:step_size:wav_size]';   %group delay is adjusted
    seg_pcm_stop{w} = min(seg_pcm_start{w}+winLenSTMSP(w),wav_size);
    seg_pcm_stop{w}(1) = min(group_delay,wav_size);
    seg_pcm_num{w} = size(seg_pcm_start{w},1);
    f_pitch_energy{w} = zeros(120,seg_pcm_num{w});
end


for p=parameter.midiMin:parameter.midiMax
%     fprintf('\b\b\b\b');fprintf('%4i',p);
    index = fs_index(p);
    f_filtfilt = filtfilt(h(p).b, h(p).a, pcm_ds{index});
    f_square = f_filtfilt.^2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % f_pitch_energy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for w=1:length(winLenSTMSP)
        factor = (parameter.fs/fs_pitch(p));      %adjustment for sampling rate
        for k=1:seg_pcm_num{w}
            start = ceil((seg_pcm_start{w}(k)/parameter.fs)*fs_pitch(p));
            stop = floor((seg_pcm_stop{w}(k)/parameter.fs)*fs_pitch(p));
            f_pitch_energy{w}(p,k)=sum(f_square(start:stop))*factor;
        end
    end
end
% fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save f_pitch_energy for each window size separately as f_pitch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sideinfo.pitch.version = 1;
sideinfo.pitch.midiMin = parameter.midiMin;
sideinfo.pitch.midiMax = parameter.midiMax;
if parameter.save == 1
    for w=1:num_window;
        f_pitch = f_pitch_energy{w};
        sideinfo.pitch.winLenSTMSP = winLenSTMSP(w);
        sideinfo.pitch.winOvSTMSP = winOvSTMSP(w);
        sideinfo.pitch.featureRate = featureRate(w);
        sideinfo.pitch.shiftFB = parameter.shiftFB;
        sideinfo.pitch.featuresAreTuned = 0;
        if parameter.saveAsTuned
            sideinfo.pitch.featuresAreTuned = 1;
            filename = strcat(parameter.saveFilename,'_pitch_',num2str(winLenSTMSP(w)));
        else
            switch(parameter.shiftFB)
                case 0
                    filename = strcat(parameter.saveFilename,'_pitch_',num2str(winLenSTMSP(w)));
                case 1
                    filename = strcat(parameter.saveFilename,'_pitch_',num2str(winLenSTMSP(w)),'_minusQuarter');
                case 2
                    filename = strcat(parameter.saveFilename,'_pitch_',num2str(winLenSTMSP(w)),'_minusThird');
                case 3
                    filename = strcat(parameter.saveFilename,'_pitch_',num2str(winLenSTMSP(w)),'_minusHalf');
                case 4
                    filename = strcat(parameter.saveFilename,'_pitch_',num2str(winLenSTMSP(w)),'_minusTwoThird');
                case 5
                    filename = strcat(parameter.saveFilename,'_pitch_',num2str(winLenSTMSP(w)),'_minusThreeQuarter');
            end
        end
        save(strcat(parameter.saveDir,filename),'f_pitch','sideinfo');
    end
else
    f_pitch = f_pitch_energy{num_window};
    sideinfo.pitch.winLenSTMSP = winLenSTMSP(num_window);
    sideinfo.pitch.winOvSTMSP = winOvSTMSP(num_window);
    sideinfo.pitch.featureRate = featureRate(num_window);
    sideinfo.pitch.shiftFB = parameter.shiftFB;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if parameter.visualize == 1
    for w=1:num_window;
        parameterVis.featureRate = featureRate(w);
        visualizePitch(f_pitch_energy{w},parameterVis);
    end
end

end


