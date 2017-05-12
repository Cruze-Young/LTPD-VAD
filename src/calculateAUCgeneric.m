function [ AUC_ALL ] = calculateAUCgeneric(results)
    nSNR = size(results, 2) - 1;
    nNoise = size(results, 1) - 1;
    nAlg = size(results{2, 2}, 1);
    AUC_ALL = cell(2, nSNR);
    for s = 1 : nSNR
        AUC = cell(nNoise+1,nAlg);

        for i = 1 : nNoise
            AUC{i+1,1} = results{i+1, 1};
        end
        AUC{i+2,1} = 'average';

        for i = 1 : nAlg
            AUC{1,i+1} = results{2, s+1}{i, 1};
        end

        auc = zeros(nNoise, nAlg);
        for i = 1 : nNoise
            auc(i,:) = singleAUC(results{i+1,s+1});
        end

        for i = 1:nNoise
            for j = 1:nAlg
                AUC{i+1,j+1} = auc(i,j);
            end
        end

        AUC = AUC';

        for i = 2:(nAlg+1)
            tmp = 0;
            for j = 2 : nNoise+1
                tmp = tmp + AUC{i,j};
            end
            AUC{i,nNoise+2} = (1/nNoise)*tmp;
        end
        
        AUC_ALL{1, s} = results{1, s + 1};
        AUC_ALL{2, s} = AUC;
    end
end

function [ AUC ] = singleAUC(results)
    nAlg = size(results, 1);
    AUC = zeros(1,nAlg);
    for i = 1 : nAlg
        AUC(i) = trapz(flipud([1;results{i,2}(:,3);0]),flipud([1;results{i,2}(:,2);0]));
    end
end