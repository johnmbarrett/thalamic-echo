figure
datas = {Response_start_2 Response_peak_2 Response_vol_2};
dataLabels = {'Response start' 'Response peak' 'Response volume'};

for ii = 1:3
    cc = [min(datas{ii}(isfinite(datas{ii}))) max(datas{ii}(isfinite(datas{ii})))];
    
    for jj = 1:4
        subplot(4,4,4*(ii-1)+jj);
        imagesc(datas{ii}(:,:,jj));
        caxis(cc);
        colormap(hot);
        
        if ii == 1
            title(sprintf('Probe %d',Probes(jj)));
        end
        
        if jj == 1
            ylabel(dataLabels{ii});
        end
    end
end

probeLocations = [10 7; 9 5; 4 5;  3 7; NaN NaN]-1;

d = sqrt(bsxfun(@minus,params.X',probeLocations(:,1)).^2+bsxfun(@minus,11-params.Y',probeLocations(:,2)).^2);
dd = zeros(12,12,5);

for ii = 1:4
    subplot(4,4,12+ii);
    dd(:,:,ii) = createMap(d(Probes(ii),:),params);
    imagesc(dd(:,:,ii));
end

%%

noResponse = bsxfun(@eq,Response_start_2,permute(max(reshape(Response_start_2,[],size(Response_start_2,3))),[3 1 2])); %Response_vol_2 == 0; %latencyIndex == 0;
% datas = {volume peak latencyIndex/1000};
% dataLabels = {'Response Volume' 'Response Peak' 'Threshold Latency'};
% probes = [1:3 5];
nProbes = numel(Probes);

beta = zeros(3,nProbes,3,2);
R2 = zeros(3,nProbes,3); % datas, probes, hemispheres

figure
for ii = 1:3
    datas{ii}(noResponse) = NaN;
    yy = [0 max(datas{ii}(isfinite(datas{ii})))];
    
    for jj = 1:nProbes
        subplot(3,nProbes,nProbes*(ii-1)+jj);
        hold on;
        
        x = cell(1,2);
        y = cell(1,2);
        
        for kk = 1:2
            x{kk} = reshape(dd(:,(6*(kk-1)+1):6*kk,Probes(jj)),[],1);
            y{kk} = reshape(datas{ii}(:,(6*(kk-1)+1):6*kk,Probes(jj)),[],1);
            plot(x{kk},y{kk},'Color',[kk-1 0 2-kk],'LineStyle','none','Marker','o');
            [beta(ii,jj,kk,:),~,~,~,stats] = regress(y{kk}(isfinite(y{kk})),[ones(sum(isfinite(y{kk})),1) x{kk}(isfinite(y{kk}))]);
            plot(0:10,beta(ii,jj,kk,1)+(0:10)*beta(ii,jj,kk,2),'Color',3*[kk-1 0 2-kk]/4);
            R2(ii,jj,kk) = stats(4);
        end
        
        x = vertcat(x{:});
        y = vertcat(y{:});
        [beta(ii,jj,3,:),~,~,~,stats] = regress(y(isfinite(y)),[ones(sum(isfinite(y)),1) x(isfinite(y))]);
        plot(0:10,beta(ii,jj,3,1)+(0:10)*beta(ii,jj,3,2),'Color',[0.75 0 0.75]);
        R2(ii,jj,3) = stats(4);
        7
        ylim(yy);
        
        if ii == 1
            title(sprintf('Probe %d',Probes(jj)));
        end
        
        if jj == 1
            ylabel(dataLabels{ii});
        end
    end
end