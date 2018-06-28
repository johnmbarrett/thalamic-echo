function plotPSTH(psth,params,isParametric)
    assert(ismember(ndims(psth),[3 4]),'PSTH should be a trial-average or grand average.');
    assert(istable(params) && all(ismember({'X' 'Y' 'PulseWidth' 'PercentPower'},params.Properties.VariableNames)),'Params should be a table with fields X, Y, PercentPower, and PulseWidth');
    
    if nargin < 3
        isParametric = false;
    end

    if isParametric
        [~,~,figureIndex] = unique([params.X params.Y],'rows');
        rowField = 'PulseWidth';
        colField = 'PercentPower';
    else
        probeIndices = params.X == 100 | params.Y == 100;
        index = [repmat({':'},1,ndims(psth)-1) ~probeIndices];
        psth = psth(index{:});
        params = params(~probeIndices,:);
        [~,~,figureIndex] = unique([params.PulseWidth params.PercentPower],'rows');
        rowField = 'Y';
        colField = 'X';
    end
    
    rowValues = unique(params.(rowField));
    colValues = unique(params.(colField));
    
    minPSTH = min(psth(:));
    maxPSTH = max(psth(:));
    
    assert(isfinite(minPSTH) && isfinite(maxPSTH),'All PSTH values must be finite');
    
    if minPSTH == maxPSTH
        warning('thalamic-echo:FlatPSTH','PSTH is flat.'); %#ok<CTPCT>
        maxPSTH = minPSTH+1;
    end
    
    rows = numel(rowValues);
    cols = numel(colValues);
    timeIndex = 50:200; % TODO : specify time limits
    t = (timeIndex-100)/1000;
    
    isGrandAverage = ndims(psth) == 3;
    
    if ~isGrandAverage
        % re-stack probes
        psth = reshape(psth,size(psth,1),[],size(psth,4));
    end
    
    for ii = 1:max(figureIndex)
        figure
        
        for jj = 1:rows
            for kk = 1:cols
                subplot(rows,cols,cols*(jj-1)+kk);
                
                paramIndex = params.(rowField) == rowValues(jj) & params.(colField) == colValues(kk) & figureIndex == ii;
                
                if isGrandAverage
                    plot(t,psth(timeIndex,:,paramIndex));
                    xlim(t([1 end]));
                    ylim([minPSTH maxPSTH]);
                else
                    imagesc(psth(timeIndex,:,paramIndex)')
                    colormap(jet);
                    caxis([minPSTH maxPSTH]);
                end
            end
        end
    end
end