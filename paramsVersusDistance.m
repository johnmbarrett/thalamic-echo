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
        [Y,X] = ndgrid(1:x,1:y);
        
        load('response_params.mat');
    
        figure
        datas = {Response_start_2 Response_peak_2 Response_vol_2};
        dataLabels = {'Response start' 'Response peak' 'Response volume'};
        probes = 1:size(Response_start_2,3);
        nProbes = numel(probes); %%min(4,numel(probes));
        nCortexProbes = min(4,nProbes);

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
        
        probeX = zeros(1,nProbes);
        probeY = zeros(1,nProbes);

        for ii = 1:nCortexProbes
            subplot(4,nProbes,(3*nProbes)+ii);
            [probeX(ii),probeY(ii)] = ind2sub([y x],probeLocations{gg}(hh,ii));
            dd(:,:,ii) = sqrt((Y-probeY(ii)).^2+(X-probeX(ii)).^2); %createMap(d(probes(ii),:),params);
            imagesc(dd(:,:,ii));
            daspect([1 1 1]);
        end
        
        m = ceil(size(dd,2)/2);
        side = ceil(probeX/m);
        
        annotation('textbox', [0 0.9 1 0.1], 'String', sprintf('Date %s recording %s',dates(gg,:),mps{gg}(hh).name), 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Interpreter', 'none')
        
        jbsavefig(gcf,'all_maps_%s_%s',dates(gg,:),mps{gg}(hh).name);
%         close(gcf);

        %%

        noResponse = bsxfun(@eq,Response_start_2,permute(max(reshape(Response_start_2,[],size(Response_start_2,3))),[3 1 2])); %Response_vol_2 == 0; %latencyIndex == 0;
        % datas = {volume peak latencyIndex/1000};
        % dataLabels = {'Response Volume' 'Response Peak' 'Threshold Latency'};
        % probes = [1:3 5];

        beta = zeros(3,nProbes,4,2);
        R2 = zeros(3,nProbes,4); % datas, probes, hemispheres

        figure
        for ii = 1:3
            datas{ii}(noResponse) = NaN;
            yy = [0 max(datas{ii}(isfinite(datas{ii})))];

            for jj = 1:nCortexProbes
                subplot(4,nProbes,nProbes*(ii-1)+jj);
                hold on;

                x = cell(1,2);
                y = cell(1,2);

                for kk = 1:2
                    x1 = 1+(kk-1)*ceil(m);
                    x2 = ceil(m)*kk-mod(size(dd,2),2);
                    sign = (3-2*kk);
                    x{kk} = sign*reshape(dd(:,x1:x2,probes(jj)),[],1);
                    y{kk} = reshape(datas{ii}(:,x1:x2,probes(jj)),[],1);
                    plot(x{kk},y{kk},'Color',[kk-1 0 2-kk],'LineStyle','none','Marker','o','MarkerSize',3*(1+(kk==side(jj))));
                    
                    [beta(ii,jj,kk,:),~,~,~,stats] = regress(y{kk}(isfinite(y{kk})),[ones(sum(isfinite(y{kk})),1) x{kk}(isfinite(y{kk}))]);
                    
                    minD = sign*min(abs(x{kk}(isfinite(y{kk}))));
                    maxD = sign*max(abs(x{kk}(isfinite(y{kk}))));
                    
                    if kk == side(jj)
                        plot([minD maxD],beta(ii,jj,kk,1)+[minD maxD]*beta(ii,jj,kk,2),'Color',3*[kk-1 0 2-kk]/4);
                    else
                        probeXdash = size(dd,2)-probeX(jj)+1; % TODO : flip about actual midline
                        dddash = sqrt((Y-probeY(jj)).^2+(X-probeXdash).^2);
                        xdash = sign*reshape(dddash(:,x1:x2),[],1);
                        
                        plot(xdash,y{kk},'Color',[0.75 0.75 0.75],'LineStyle','none','Marker','o','MarkerSize',3);
                        
                        [beta(ii,jj,4,:),~,~,~,stats] = regress(y{kk}(isfinite(y{kk})),[ones(sum(isfinite(y{kk})),1) xdash(isfinite(y{kk}))]);
                        
                        maxDdash = sign*max(abs(xdash(isfinite(y{kk}))));
                        
                        plot([0 maxDdash],beta(ii,jj,4,1)+[0 maxDdash]*beta(ii,jj,4,2),'Color',[0.5 0.5 0.5]);
                        
                        R2(ii,jj,4) = stats(4);
                    end
                    
                    R2(ii,jj,kk) = stats(4);
                end

                x = vertcat(x{:});
                y = vertcat(y{:});
                [beta(ii,jj,3,:),~,~,~,stats] = regress(y(isfinite(y)),[ones(sum(isfinite(y)),1) x(isfinite(y))]);
%                 plot(0:d,beta(ii,jj,3,1)+(0:d)*beta(ii,jj,3,2),'Color',[0.75 0 0.75]);
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
        
        for ii = 1:nProbes
            subplot(4,nProbes,3*nProbes+ii);
            
            hold on;
            
            l = datas{1}(:,:,probes(ii));
            v = datas{2}(:,:,probes(ii));
            plot(l,v,'Color','k','LineStyle','none','Marker','o');
            
            b = regress(v(:),[ones(numel(l),1) l(:)]);
            
            plot([0 max(l(:))],b(1)+[0 max(l(:))]*b(2),'Color',[0.5 0.5 0.5],'LineStyle','--');
            
            xlim([0 max(datas{1}(:))]);
            ylim([0 max(datas{2}(:))]);
            
            xlabel('Response latency');
            
            if ii == 1
                ylabel('Response volume');
            end
        end
        
        slopes{gg}{hh} = beta(:,:,:,2);
        intercepts{gg}{hh} = beta(:,:,:,2);
        R2s{gg}{hh} = R2;
        
        annotation('textbox', [0 0.9 1 0.1], 'String', sprintf('Date %s recording %s',dates(gg,:),mps{gg}(hh).name), 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Interpreter', 'none')
        
        jbsavefig(gcf,'responses_versus_distance_%s_%s',dates(gg,:),mps{gg}(hh).name);
%         close(gcf);
    end
end

cd(topDir);

comment = sprintf('Slope is in ms/mm\n\nOuter cell array is date; inner cell array is multiple recordings on the same day; rows are latency, then peak, then vol; columns are probes; pages are left hemisphere, then right, then bilateral');

save('spread_parameters','comment','dates','mps','slopes','intercepts','R2s');