clc; clear; close all;
addpath ./src
addpath ./src/voicebox;
addpath ./src/chroma
%%
load ./data/DRx.mat;
load ./data/noiseAll.mat
load ./data/thr.mat
vadinfo(1, :) = [{@LTSDVAD}, {thrLTSD}];
vadinfo(2, :) = [{@LTSVVAD}, {thrLTSV}];
vadinfo(3, :) = [{@LTPDVAD}, {thrLTPD}];
vadinfo(4, :) = [{@sohnVAD}, {thrSohn}];
vadinfo(5, :) = [{@harmfreqVAD}, {thrHarm}];
vadfns = {'LTSDVAD'; 'LTSVVAD'; 'LTPDVAD'; 'sohnVAD'; 'harmfreqVAD'};

Nx = {babble; buccaneer1; buccaneer2; destroyerengine; destroyerops;...
      f16; factory1; factory2; hfchannel; m109;...
      leopard; machinegun; pink; volvo; white};
%%
fs = 16000;
wsec = 0.05;
isec = 0.05;
SNR = [-5, 0, 5, 10, 15, 20];
results = [{''}; {'babble'}; {'buccaneer1'}; {'buccaneer2'}; {'destroyerengine'}; ...
           {'destroyerops'}; {'f16'}; {'factory1'}; {'factory2'}; {'hfchannel'};...
           {'m109'}; {'leopard'}; {'machinegun'}; {'pink'}; {'volvo'}; {'white'}]; 
results{1,2} = '-5dB'; 
results{1,3} = '0dB'; 
results{1,4} = '5dB'; 
results{1,5} = '10dB'; 
results{1,6} = '15dB'; 
results{1,7} = '20dB'; 

parameter.N = 3;
parameter.R = 6;
parameter.M = 2;
for i = 2 : size(results, 1)
    display(strcat('Noise (',num2str(i-1),'/',num2str(length(Nx)),'): ',results{i,1}));
    for j = 2 : size(results, 2)
        display(strcat('SNR:',num2str(SNR(j-1)), 'dB'));
        results{i,j} = evaluateROC(DRx, Nx{i-1}, fs, SNR(j-1), wsec, isec, vadinfo, parameter);
    end
end
