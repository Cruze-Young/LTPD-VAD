function plotROCgeneric(results, vadfns)
    nNoise = size(results, 1) - 1;
    nRow = ceil(nNoise / 2);
    if nNoise > 1
        nCol = 2;
    else
        nCol = 1;
    end
    nfigure = size(results, 2) - 1;
    for f = 2 : nfigure+1
        figure(f);
        for r = 2 : nNoise+1
            subplot(nRow,nCol,r-1); 
            plotLine(results{r, f}); 
            title(strcat(results{r,1},{' '},results{1,f})); 
            ylabel('True Positive Rate'); 
            xlabel('False Positive Rate');
        end
        legend(vadfns);
    end
end

function plotLine(results)
    hold on;
    nLine = size(results, 1);
    Color = [0 0 1; 1 0 0; 1 0 1; ... 
                 0 1 1; 0 1 0; 0 0 0; ...
                 0.4 0.7 0.4; 0.7 0.4 0.2; 0.8 0.3 0.6];
    lineStyle = { '-+'; '-o'; '-*'; '-.'; '-x'; '-square'; '-diamond'; '-v'; '-^'; ...
                  '->';  '-<'; '-pentagram'; '-hexagram' };
    for l = 1 : nLine
        plot([1;results{l,2}(:,3);0],[1;results{l,2}(:,2);0], lineStyle{l}, 'Color', Color(l, :));
    end
end
