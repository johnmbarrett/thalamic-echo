function [volume,peak,peakIndex,latencyIndex] = calculateResponseParameters(mpsth,params,psthFigure,isPlot)
    stimulusIndex = 101; % TODO : specify?
    baseline = mpsth(1:(stimulusIndex-1),:,:);
    % TODO : per probe or even per channel baseline is better but 
    % eyeballing the plots it doesn't make too much difference
    mu = mean(baseline(:));
    sigma = std(baseline(:));
    
    rows = numel(unique(params.Y));
    cols = numel(unique(params.X));
    
    if nargin > 2 && isgraphics(psthFigure)
        styles = {'-' '--' ':'};
        for jj = 1:rows
            for kk = 1:cols
                subplot(rows,cols,cols*(jj-1)+kk);
                
                hold on;
                
                for ll = 1:3
                    line(repmat(xlim',1,2),repmat(mu+[-1 1]*(ll+3)*sigma,2,1),'Color','k','LineStyle',styles{ll});
                end
            end
        end     
    end
    
    responseLength = 100; % TODO : specify
    responses = mpsth(stimulusIndex:(stimulusIndex+responseLength-1),:,:);
    
    volume = squeeze(sum(responses));
    
    [peak,peakIndex] = max(responses);
    peak = squeeze(peak);
    peakIndex = squeeze(peakIndex);
    
    latencyIndex = squeeze(cellfun(@(r) sum(find(r > mu+3*sigma,1)),mat2cell(responses,size(responses,1),ones(1,size(responses,2)),ones(1,size(responses,3)))));
    
    if (nargin > 4 && ~isPlot) || (nargin > 3 && islogical(psthFigure) && ~psthFigure)
        return % plot by default
    end
    
    figure
    datas = {volume peak peakIndex/1000 latencyIndex/1000};
    dataLabels = {'Response Volume' 'Response Peak' 'Peak Latency' 'Threshold Latency'};
    nDatas = numel(datas);
    nProbes = size(volume,1);
    
    for ii = 1:nDatas
        clim = [min(datas{ii}(isfinite(datas{ii}))) max(datas{ii}(isfinite(datas{ii})))];
        
        if clim(1) == clim(2)
            clim(2) = clim(1)+1;
        end
        
        for jj = 1:nProbes
            subplot(nDatas,nProbes,nProbes*(ii-1)+jj);
            
            imagesc(createMap(datas{ii}(jj,:),params));
            
            caxis(clim);
            colorbar;
            colormap(hot);
            
            if ii == 1
                title(sprintf('Probe %d',jj));
            end
            
            if jj == 1
                ylabel(dataLabels{ii});
            end
        end
    end
end