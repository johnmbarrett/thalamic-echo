% TODO : this preamble should probably be refactored into a function/script

topDir = 'Z:\LACIE\Manuscripts\2018 in vivo LSPS Ntsr1 etc\data';
cd(topDir);

%load('probe_locations.mat');
strain = 'SepW'; %'Ntsr1';
experimentParams = readExperimentSpreadsheet(strain);
% experimentParams.Midline = cellfun(@str2num,experimentParams.Midline,'UniformOutput',false);

if ~strcmp(strain,'Ntsr1')
    topDir = [topDir '\' strain];
    cd(topDir)
end

%%

nRecordings = size(experimentParams,1);
slopes = cell(nRecordings,1);
intercepts = cell(nRecordings,1);
R2s = cell(nRecordings,1);

colormaps = {jet2(256) hot(256) hot(256)};
colormaps{1}(2:end,:) = flipud(colormaps{1}(2:end,:));

probeNames = {'R-S1' 'R-M1' 'L-M1' 'L-S1' 'Thal'};

edges = 0:1:10;
binCentres = mean([edges(1:end-1);edges(2:end)]);
nBins = numel(edges)-1;

binnedData = cell(nBins,2,4,3,nRecordings);

for gg = 1:nRecordings
    cd([topDir '\' datestr(experimentParams.Date(gg),'yyyymmdd') '\' experimentParams.MPFolder{gg}]);

    x = experimentParams.X(gg);
    y = experimentParams.Y(gg);
    [Y,X] = ndgrid(1:x,1:y);

    load('response_params.mat');

    figure
    datas = {Response_start_2 Response_peak_2 Response_vol_2};
    dataLabels = {'Response start' 'Response peak' 'Response volume'};
    probes = 1:size(Response_start_2,3);
    nProbes = numel(probes); %%min(4,numel(probes));
    nCortexProbes = min(4,nProbes);

    noResponse = bsxfun(@eq,Response_start_2,permute(max(reshape(Response_start_2,[],size(Response_start_2,3))),[3 1 2])); %Response_vol_2 == 0; %latencyIndex == 0;
    
    probeLocations = arrayfun(@(ii) experimentParams(gg,:).(sprintf('ProbeLocations_%d',ii)),1:4,'UniformOutput',false);
    probeLocations = [probeLocations{:}];
    
    % need these for indexing into the name array
    validProbeIndices = find(isfinite(probeLocations));
    
    assert(numel(validProbeIndices) == nCortexProbes);
    
    probeLocations = probeLocations(validProbeIndices);

    for ii = 1:3
        datas{ii}(noResponse) = NaN;
        datas{ii}(experimentParams.ExcludeRows{gg},:,:) = NaN;
        datas{ii}(:,experimentParams.ExcludeColumns{gg},:) = NaN;
        
        m = experimentParams.Midline{gg};
        
        if ~isscalar(m) || round(m) == m
            datas{ii}(:,m,:) = NaN;
        end
        
        for jj = 1:size(experimentParams.ExcludePixels{gg},1)
            datas{ii}(experimentParams.ExcludePixels{gg}(jj,1),experimentParams.ExcludePixels{gg}(jj,2),:) = NaN;
        end
        
        cc = [min(datas{ii}(isfinite(datas{ii}))) max(datas{ii}(isfinite(datas{ii})))];

        for jj = 1:nProbes
            subplot(4,nProbes,nProbes*(ii-1)+jj);
            imagesc(datas{ii}(:,:,jj));
            caxis(cc);
            colormap(gca,colormaps{ii});

            if ii == 1
                if jj <= nCortexProbes
                    title(probeNames{validProbeIndices(jj)});
                else
                    title(probeNames{end});
                end
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
        [probeX(ii),probeY(ii)] = ind2sub([y x],probeLocations(ii));
        dd(:,:,ii) = sqrt((Y-probeY(ii)).^2+(X-probeX(ii)).^2); %createMap(d(probes(ii),:),params);
        imagesc(dd(:,:,ii));
        daspect([1 1 1]);
        colormap(gca,gray(256));
    end

    side = (probeX > m(end))+1; % TODO : hopefully none of the probes are in the midline 

    if nProbes > nCortexProbes
        side(end) = 2; % TODO : always R-thalamus?
    end

    annotation('textbox', [0 0.9 1 0.1], 'String', sprintf('Date %s recording %s',datestr(experimentParams.Date(gg),'yyyymmdd'),experimentParams.MPFolder{gg}), 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Interpreter', 'none')

    jbsavefig(gcf,'%s\\plots\\all_maps_%s_%s',topDir,datestr(experimentParams.Date(gg),'yyyymmdd'),experimentParams.MPFolder{gg});
%         close(gcf);

    %%
    % datas = {volume peak latencyIndex/1000};
    % dataLabels = {'Response Volume' 'Response Peak' 'Threshold Latency'};
    % probes = [1:3 5];

    beta = zeros(3,nProbes,4,2);
    R2 = zeros(3,nProbes,4); % datas, probes, hemispheres
    [left,right] = getHemispheres(size(dd,2),m);
    hemispheres = {left right};

    figure
    for ii = 1:3
        yy = [0 max(datas{ii}(isfinite(datas{ii})))];

        for jj = 1:nCortexProbes
            subplot(4,nProbes,nProbes*(ii-1)+jj);
            hold on;

            x = cell(1,2);
            y = cell(1,2);

            for kk = 1:2
                sign = 2*kk-3;
                x{kk} = sign*reshape(dd(:,hemispheres{kk},probes(jj)),[],1);
                y{kk} = reshape(datas{ii}(:,hemispheres{kk},probes(jj)),[],1);
                plot(x{kk},y{kk},'Color',[kk-1 0 2-kk],'LineStyle','none','Marker','o','MarkerSize',3*(1+(kk==side(jj))));

                [beta(ii,jj,kk,:),~,~,~,stats] = regress(y{kk}(isfinite(y{kk})),[ones(sum(isfinite(y{kk})),1) x{kk}(isfinite(y{kk}))]);

                minD = sign*min(abs(x{kk}(isfinite(y{kk}))));
                maxD = sign*max(abs(x{kk}(isfinite(y{kk}))));

                if kk == side(jj)
                    plot([minD maxD],beta(ii,jj,kk,1)+[minD maxD]*beta(ii,jj,kk,2),'Color',3*[kk-1 0 2-kk]/4);
                    
                    binnedData(:,1,validProbeIndices(jj),ii,gg) = arrayfun(@(l,r) y{kk}(abs(x{kk}) >= l & abs(x{kk}) < r),edges(1:end-1),edges(2:end),'UniformOutput',false); % TODO : last edge?
                else
                    probeXdash = probeX(jj)-2*(probeX(jj)-median(m));
                    dddash = sqrt((Y-probeY(jj)).^2+(X-probeXdash).^2);
                    xdash = sign*reshape(dddash(:,hemispheres{kk}),[],1);

                    plot(xdash,y{kk},'Color',[0.75 0.75 0.75],'LineStyle','none','Marker','o','MarkerSize',3);
                    
                    binnedData(:,2,validProbeIndices(jj),ii,gg) = arrayfun(@(l,r) y{kk}(abs(xdash) >= l & abs(xdash) < r),edges(1:end-1),edges(2:end),'UniformOutput',false);

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
                title(probeNames{validProbeIndices(jj)});
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

    slopes{gg} = beta(:,:,:,2);
    intercepts{gg} = beta(:,:,:,2);
    R2s{gg} = R2;

    annotation('textbox', [0 0.9 1 0.1], 'String', sprintf('Date %s recording %s',datestr(experimentParams.Date(gg),'yyyymmdd'),experimentParams.MPFolder{gg}), 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Interpreter', 'none')

    jbsavefig(gcf,'%s\\plots\\responses_versus_distance_%s_%s',topDir,datestr(experimentParams.Date(gg),'yyyymmdd'),experimentParams.MPFolder{gg});
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
                        xrange = left;
                        yrange = right;
                    else
                        yrange = hemispheres{side(ll)};
                        xrange = hemispheres{3-side(ll)};
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

    annotation('textbox', [0 0.9 1 0.1], 'String', sprintf('Date %s recording %s',datestr(experimentParams.Date(gg),'yyyymmdd'),experimentParams.MPFolder{gg}), 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Interpreter', 'none')

    jbsavefig(gcf,'%s\\plots\\laterality_%s_%s',topDir,datestr(experimentParams.Date(gg),'yyyymmdd'),experimentParams.MPFolder{gg});
    %%
    figure
    for ii = 1:3
        for jj = 1:nCortexProbes
            subplot(3,nCortexProbes,nCortexProbes*(ii-1)+jj);
            hold on;
            
            x = vertcat(binnedData{:,:,jj,ii,gg});
            
            if side(jj) == 1
                allBinCentres = [-binCentres';binCentres'];
            else
                allBinCentres = [binCentres';-binCentres'];
            end
            
            g = arrayfun(@(kk,data) repmat(kk,numel(data{1}),1),allBinCentres,reshape(binnedData(:,:,jj,ii,gg),[],1),'UniformOutput',false);
            
            boxplot(x,vertcat(g{:}),'Positions',sort(allBinCentres(~cellfun(@isempty,g))));
            
            xtick = [fliplr(-edges(2:end)) edges];
            set(gca,'XTick',xtick,'XTickLabel',xtick);
            
            xlim([-edges(end) edges(end)]);
            
            for kk = 1:2
                for ll = 1:nBins
                    if isempty(g{nBins*(kk-1)+ll})
                        continue
                    end
                    
                    plot((2*side(jj)-3)*(3-2*kk)*binCentres(ll),binnedData{ll,kk,jj,ii,gg},'Color',[0.5 0.5 0.5]*(kk==2),'LineStyle','none','Marker','o','MarkerSize',3+3*(kk==1));
                end
            end
            
            if ii == 1
                title(probeNames{validProbeIndices(jj)});
            end

            if jj == 1
                if ii == 3
                    xlabel('Distance from midline (laser grid steps)');
                end
                
                ylabel(dataLabels{ii});
            end
            
            yy = [min(datas{ii}(:)) max(datas{ii}(:))];
            yy = yy+([-1 1]/10)*diff(yy);
            
            ylim(yy);
        end
    end
    
    annotation('textbox', [0 0.9 1 0.1], 'String', sprintf('Date %s recording %s',datestr(experimentParams.Date(gg),'yyyymmdd'),experimentParams.MPFolder{gg}), 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Interpreter', 'none')

    jbsavefig(gcf,'%s\\plots\\binned_data_%s_%s',topDir,datestr(experimentParams.Date(gg),'yyyymmdd'),experimentParams.MPFolder{gg});
%%        
end
%%
close all

cd(topDir);

comment = sprintf('Slope is in ms/mm\n\nFor slopes/intercepts/R2s, outer cell array is date; inner cell array is multiple recordings on the same day; rows are latency, then peak, then vol; columns are probes; pages are left hemisphere, then right, then bilateral\n\nFor binnedData, first dimension is bins, second is ipsi/contra, third is probes, fourth is start/peak/vol, fifth is recordings.');

save('spread_parameters','comment','slopes','intercepts','R2s','binnedData','edges','binCentres','nBins');