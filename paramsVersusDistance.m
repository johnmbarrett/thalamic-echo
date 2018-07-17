% TODO : laterality corrected distance, gordon's plots

topDir = 'Z:\LACIE\Manuscripts\2018 in vivo LSPS Ntsr1 etc\data';
cd(topDir);

load('probe_locations.mat');

%%

slopes = cellfun(@(A) cell(size(A)),mps,'UniformOutput',false);
intercepts = cellfun(@(A) cell(size(A)),mps,'UniformOutput',false);
R2s = cellfun(@(A) cell(size(A)),mps,'UniformOutput',false);

for gg = 1:numel(probeLocations)
    for hh = 1:numel(mps{gg})
        cd([topDir '\' dates(gg,:) '\' mps{gg}(hh).name]);
        
        x = xy{gg}(hh,1);
        y = xy{gg}(hh,2);
        [X,Y] = ndgrid(1:x,1:y);
        
        load('response_params.mat');
    
        figure
        datas = {Response_start_2 Response_peak_2 Response_vol_2};
        dataLabels = {'Response start' 'Response peak' 'Response volume'};
        probes = 1:size(Response_start_2,3);
        nProbes = min(4,numel(probes));

        for ii = 1:3
            cc = [min(datas{ii}(isfinite(datas{ii}))) max(datas{ii}(isfinite(datas{ii})))];

            for jj = 1:nProbes
                subplot(4,nProbes,nProbes*(ii-1)+jj);
                imagesc(datas{ii}(:,:,jj));
                caxis(cc);
                colormap(hot);

                if ii == 1
                    title(sprintf('Probe %d',probes(jj)));
                end

                if jj == 1
                    ylabel(dataLabels{ii});
                end
                
                daspect([1 1 1]);
            end
        end

%         d = sqrt(bsxfun(@minus,params.X',probeLocations(:,1)).^2+bsxfun(@minus,11-params.Y',probeLocations(:,2)).^2);
        dd = zeros(x,y,nProbes);

        for ii = 1:nProbes
            subplot(4,nProbes,(3*nProbes)+ii);
            [probeY,probeX] = ind2sub([y x],probeLocations{gg}(hh,ii));
            dd(:,:,ii) = sqrt((X-probeX).^2+(Y-probeY).^2); %createMap(d(probes(ii),:),params);
            imagesc(dd(:,:,ii));
            daspect([1 1 1]);
        end
        
        d = max(dd(:));
        m = ceil(size(dd,2)/2);
        
        jbsavefig(gcf,'all_maps_%s_%s',dates(gg,:),mps{gg}(hh).name);
        close(gcf);

        %%

        noResponse = bsxfun(@eq,Response_start_2,permute(max(reshape(Response_start_2,[],size(Response_start_2,3))),[3 1 2])); %Response_vol_2 == 0; %latencyIndex == 0;
        % datas = {volume peak latencyIndex/1000};
        % dataLabels = {'Response Volume' 'Response Peak' 'Threshold Latency'};
        % probes = [1:3 5];

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
                    x1 = 1+(kk-1)*ceil(m);
                    x2 = ceil(m)*kk-mod(size(dd,2),2);
                    x{kk} = reshape(dd(:,x1:x2,probes(jj)),[],1);
                    y{kk} = reshape(datas{ii}(:,x1:x2,probes(jj)),[],1);
                    plot(x{kk},y{kk},'Color',[kk-1 0 2-kk],'LineStyle','none','Marker','o');
                    [beta(ii,jj,kk,:),~,~,~,stats] = regress(y{kk}(isfinite(y{kk})),[ones(sum(isfinite(y{kk})),1) x{kk}(isfinite(y{kk}))]);
                    plot(0:d,beta(ii,jj,kk,1)+(0:d)*beta(ii,jj,kk,2),'Color',3*[kk-1 0 2-kk]/4);
                    R2(ii,jj,kk) = stats(4);
                end

                x = vertcat(x{:});
                y = vertcat(y{:});
                [beta(ii,jj,3,:),~,~,~,stats] = regress(y(isfinite(y)),[ones(sum(isfinite(y)),1) x(isfinite(y))]);
                plot(0:d,beta(ii,jj,3,1)+(0:d)*beta(ii,jj,3,2),'Color',[0.75 0 0.75]);
                R2(ii,jj,3) = stats(4);
                
                ylim(yy);

                if ii == 1
                    title(sprintf('Probe %d',probes(jj)));
                end

                if jj == 1
                    ylabel(dataLabels{ii});
                end
            end
        end
        
        slopes{gg}{hh} = beta(:,:,:,2);
        intercepts{gg}{hh} = beta(:,:,:,2);
        R2s{gg}{hh} = R2;
        
        jbsavefig(gcf,'responses_versus_distance_%s_%s',dates(gg,:),mps{gg}(hh).name);
        close(gcf);
    end
end

cd(topDir);

comment = sprintf('Slope is in ms/mm\n\nOuter cell array is date; inner cell array is multiple recordings on the same day; rows are latency, then peak, then vol; columns are probes; pages are left hemisphere, then right, then bilateral');

save('spread_parameters','comment','dates','mps','slopes','intercepts','R2s');