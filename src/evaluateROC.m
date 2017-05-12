function [ vadinfo ] = evaluateROC(DRx,n,fs,snrr,wsec,isec,vadinfo, parameter)
    % add noise to clean speech
    ns = cell(length(DRx),1);
    for k = 1:length(DRx)
        % calculate the power per sample of speech when it is active
        speechPresent = logical(DRx{k}(:,2));
        psig = sum(DRx{k}(speechPresent,1).^2)/length(DRx{k}(speechPresent,1));
        % truncate the noise and calculate its power per sample
        imax = length(n) - size(DRx{k},1);
        istart = randi(imax, 1);
        noi = n(istart:(istart+size(DRx{k},1)-1));
        pnoi = sum(noi.^2)/length(noi);
        % calculate the scaling constant for noise
        sc = sqrt(psig/(pnoi*10^(snrr/10)));
        % scale the noise to the desired SNR
        noi = noi.*sc;
        % add noise to clean speech
        ns{k} = DRx{k}(:,1) + noi;
    end

    vadno = size(vadinfo,1);
    DRxno = length(DRx);
    % for each vad
    for i = 1:vadno
        vadfn = vadinfo{i,1};
        display(strcat('VAD ', '(', num2str(i), '/', num2str(vadno), ')'));
        thr = vadinfo{i,2};
        thrno = length(vadinfo{i,2});
        TP = zeros(1,thrno);
        TN = zeros(1,thrno);
        FP = zeros(1,thrno);
        FN = zeros(1,thrno);
        for k = 1:DRxno
            display(strcat('  DRx ', '(', num2str(k), '/', num2str(DRxno), ')'));
            % run VAD
            vad = vadfn(ns{k},fs,wsec,isec,parameter,0,thr);
            % evaluate VAD
            for j = 1 : thrno
                stats = recallPrecision(vad(j,:),DRx{k}(:,2));
                TP(j) = TP(j) + stats{1,2};
                TN(j) = TN(j) + stats{2,2};
                FP(j) = FP(j) + stats{3,2};
                FN(j) = FN(j) + stats{4,2};
            end
        end
        for j = 1 : thrno
            % true positive rate
            vadinfo{i,2}(j,2) = TP(j)/(TP(j)+FN(j));
            % false positive rate
            vadinfo{i,2}(j,3) = FP(j)/(FP(j)+TN(j));
            
            vadinfo{i,2}(j,4) = TP(j);
            vadinfo{i,2}(j,5) = TN(j);
            vadinfo{i,2}(j,6) = FP(j);
            vadinfo{i,2}(j,7) = FN(j);
        end
    end
end