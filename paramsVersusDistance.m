% TODO : laterality corrected distance, gordon's plots

topDir = 'Z:\LACIE\Manuscripts\2018 in vivo LSPS Ntsr1 etc\data';
cd(topDir);

load('probe_locations.mat');

%%

slopes = cellfun(@(A) cell(size(A)),mps,'UniformOutput',false);
intercepts = cellfun(@(A) cell(size(A)),mps,'UniformOutput',false);
R2s = cellfun(@(A) cell(size(A)),mps,'UniformOutput',false);

colormaps = {jet2(256) hot(256) hot(256)};
colormaps{1}(2:end,:) = flipud(colormaps{1}(2:end,:));

probeNames = {'R-S1' 'R-M1' 'L-M1' 'L-S1' 'Thal'};

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
        
        noResponse = bsxfun(@eq,Response_start_2,permute(max(reshape(Response_start_2,[],size(Response_start_2,3))),[3 1 2])); %Response_vol_2 == 0; %latencyIndex == 0;

        for ii = 1:3
            datas{ii}(noResponse) = NaN;
            cc = [min(datas{ii}(isfinite(datas{ii}))) max(datas{ii}(isfinite(datas{ii})))];

            for jj = 1:nProbes
                subplot(4,nProbes,nProbes*(ii-1)+jj);
                imagesc(datas{ii}(:,:,jj));
                caxis(cc);
                colormap(gca,colormaps{ii});

                if ii == 1
                    title(probeNames{jj});
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
            colormap(gca,gray(256));
        end
        
        m = size(dd,2)/2;
        side = ceil(probeX/m);
        
        if nProbes > nCortexProbes
            side(end) = 2; % TODO : always R-thalamus
        end
        
        annotation('textbox', [0 0.9 1 0.1], 'String', sprintf('Date %s recording %s',dates(gg,:),mps{gg}(hh).name), 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Interpreter', 'none')
        
        jbsavefig(gcf,'%s\\plots\\all_maps_%s_%s',topDir,dates(gg,:),mps{gg}(hh).name);
%         close(gcf);

        %%
        % datas = {volume peak latencyIndex/1000};
        % dataLabels = {'Response Volume' 'Response Peak' 'Threshold Latency'};
        % probes = [1:3 5];

        beta = zeros(3,nProbes,4,2);
        R2 = zeros(3,nProbes,4); % datas, probes, hemispheres

        figure
        for ii = 1:3
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
                    title(probeNames{jj});
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
            v = datas{3}(:,:,probes(ii));
            plot(l,v,'Color','k','LineStyle','none','Marker','o');
            
            b = regress(v(:),[ones(numel(l),1) l(:)]);
            
            plot([min(l(:)) max(l(:))],b(1)+[min(l(:)) max(l(:))]*b(2),'Color',[0.5 0.5 0.5],'LineStyle','--');
            
            xlim([0 max(datas{1}(:))]);
            ylim([0 max(datas{3}(:))]);
            
            xlabel('Response latency');
            
            if ii == 1
                ylabel('Response volume');
            end
        end
        
        slopes{gg}{hh} = beta(:,:,:,2);
        intercepts{gg}{hh} = beta(:,:,:,2);
        R2s{gg}{hh} = R2;
        
        annotation('textbox', [0 0.9 1 0.1], 'String', sprintf('Date %s recording %s',dates(gg,:),mps{gg}(hh).name), 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Interpreter', 'none')
        
        jbsavefig(gcf,'%s\\plots\\responses_versus_distance_%s_%s',topDir,dates(gg,:),mps{gg}(hh).name);
%         close(gcf);
%%
        figure
        colours = 'rgbmc';
        dataNames = {'amplitude' 'latency'};
        labels = {'Left ' 'Right '; 'Contralateral ' 'Ipsilateral '};
        
        for ii = 1:2
            for jj = 1:2
                for kk = 1:2
                    subplot(2,4,ii+2*(jj-1)+4*(kk-1));
                    hold on;

                    for ll = 1:max(1,(ii==2)*nProbes)
                        if kk == 1
                            xrange = 1:floor(m);
                            yrange = (ceil(m)+1):size(dd,2);
                        else
                            yrange = ((1+ceil(m)*(side(ll)-1))):ceil(m)*side(ll)-mod(size(dd,2),2);
                            xrange = ((1+ceil(m)*(2-side(ll)))):ceil(m)*(3-side(ll))-mod(size(dd,2),2);
                        end
                        
                        plot(datas{5-2*jj}(:,xrange,ll),datas{5-2*jj}(:,yrange,ll),'Color',colours(ll),'LineStyle','none','Marker','o','MarkerSize',3);
                    end
                    
                    if prod([ii jj kk] == 1)
                        hs = gobjects(nProbes,1);
                        
                        for ll = 1:nProbes
                            hs(ll) = plot(NaN,NaN,'Color',colours(ll),'LineStyle','none','Marker','o','MarkerSize',3);
                        end
                        
                        legend(hs,probeNames,'AutoUpdate','off','Location','Best');
                    end
                        
                    plot([1e-3 1e3],[1e-3 1e3],'Color',[0.5 0.5 0.5],'LineStyle','--');
                    
                    if kk == 1
                        set(gca,'XScale','log','YScale','log');
                    end
                    
                    pbaspect([1 1 1]);
                    
                    xlabel([labels{kk,1} dataNames{jj}]);
                    xlim([(kk==1)*min(datas{5-2*jj}(:)) max(datas{5-2*jj}(:))]);
                    ylabel([labels{kk,2} dataNames{jj}]);
                    ylim([(kk==1)*min(datas{5-2*jj}(:)) max(datas{5-2*jj}(:))]);
                end
            end
        end
        
        annotation('textbox', [0 0.9 1 0.1], 'String', sprintf('Date %s recording %s',dates(gg,:),mps{gg}(hh).name), 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Interpreter', 'none')
        
        jbsavefig(gcf,'%s\\plots\\laterality_%s_%s',topDir,dates(gg,:),mps{gg}(hh).name);
%%        
    end
end

close all

cd(topDir);

comment = sprintf('Slope is in ms/mm\n\nOuter cell array is date; inner cell array is multiple recordings on the same day; rows are latency, then peak, then vol; columns are probes; pages are left hemisphere, then right, then bilateral');

save('spread_parameters','comment','dates','mps','slopes','intercepts','R2s');