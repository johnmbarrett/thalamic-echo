topDir = 'Z:\LACIE\Manuscripts\2018 in vivo LSPS Ntsr1 etc\data';
cd(topDir);

load('spread_parameters.mat');

probeNames = {'R-S1' 'R-M1' 'L-M1' 'L-S1' 'Thal'};
dataLabels = {'Response start' 'Response peak' 'Response volume'};

%%

medianData = permute(cellfun(@nanmedian,binnedData),[5 1:4]); % bring the reps to the first dimension
nRecordings = size(medianData,1);
n = permute(cellfun(@(A) sum(isfinite(A)),binnedData),[5 1:4]);
N = squeeze(sum(n > 0,1));

assert(isequal(isfinite(medianData),n ~= 0));

side = [2 2 1 1]; % TODO : don't hard code this

%%

plotOptions = {'lines'};
% plotOptions = {'boxplot' 'lines'};

% TODO : there's a lot of commonality between all these different plot
% loops, that could probably be refactored out but I'm too lazy and it's
% low priority right now
for hh = 1:numel(plotOptions)
    figure
    for ii = 1:3
        for jj = 1:4
            subplot(3,4,4*(ii-1)+jj);
            hold on;

            switch plotOptions{hh}
                case 'boxplot'
                    x = reshape(medianData(:,:,:,jj,ii),[],1);

                    if side(jj) == 1
                        allBinCentres = [-binCentres';binCentres'];
                    else
                        allBinCentres = [binCentres';-binCentres'];
                    end

                    g = kron(allBinCentres,ones(6,1));

                    boxplot(x,g,'Positions',unique(g));

                    xtick = [fliplr(-edges(2:end)) edges];
                    set(gca,'XTick',xtick,'XTickLabel',xtick);

                    for kk = 1:2
                        for ll = 1:nBins
                            plot((2*side(jj)-3)*(3-2*kk)*binCentres(ll),medianData(:,ll,kk,jj,ii),'Color',[0.5 0.5 0.5]*(kk==2),'LineStyle','none','Marker','o','MarkerSize',3+3*(kk==1));
                        end
                    end
                case 'lines'
                    for kk = 1:2
                        set(gca,'ColorOrderIndex',1);

                        hs = plot((2*side(jj)-3)*(3-2*kk)*binCentres,medianData(:,:,kk,jj,ii),'Marker','o','MarkerSize',3+3*(kk==1));
                        plot((2*side(jj)-3)*(3-2*kk)*binCentres,nanmedian(medianData(:,:,kk,jj,ii)),'Marker','o','MarkerSize',3+3*(kk==1),'Color',(kk==2)*[1 1 1]/2,'LineWidth',2);

                        for ll = 1:numel(hs)
                            set(hs(ll),'Color',(1-kk/4)*get(hs(ll),'Color')+(kk==2)/2);
                        end
                    end
                otherwise
                    error('Unknown plot option.');
            end

            xlim([-edges(end) edges(end)]);

            if ii == 1
                title(probeNames{jj});
            end

            if jj == 1
                if ii == 3
                    xlabel('Distance from midline (laser grid steps)');
                end

                ylabel(dataLabels{ii});
            end

            ylim([0 max(reshape(medianData(:,:,:,:,ii),[],1))]);
        end
    end
    
    annotation('textbox',[0 0.9 1 0.1],'String',sprintf('Group data - %s',plotOptions{hh}),'EdgeColor','none','HorizontalAlignment','center','Interpreter','none')
    
    jbsavefig(gcf,'group_data_%s',plotOptions{hh});
end

%%

lmes = cell(1,3);
glmes = cell(1,3);

nBins = size(binnedData,1);
distance = repmat((1:nBins)',1,2,4,nRecordings);
laterality = repmat('ic',nBins,1,4,nRecordings);
side = repmat(reshape('rrll',1,1,4),nBins,2,1,nRecordings);
region = repmat(reshape('smms',1,1,4),nBins,2,1,nRecordings);
mouse = repmat(reshape(1:nRecordings,1,1,1,nRecordings),nBins,2,4,1);

groups = {distance laterality side region mouse};

dvNames = {'Latency' 'Peak' 'Volume'};
ivNames = {'Distance' 'Laterality' 'Side' 'Region' 'Mouse'};

distributions = {'Gamma' 'Poisson' 'Poisson'};

for ii = 1:3
    dv = cellfun(@(A) A(isfinite(A))',binnedData(:,:,:,1,:),'UniformOutput',false);
    dv = [dv{:}]';
    
    ivs = cell(1,5);
    
    for jj = 1:5
        ivs{jj} = arrayfun(@(A,B) repmat(A,sum(isfinite(B{1})),1),groups{jj},squeeze(binnedData(:,:,:,ii,:)),'UniformOutput',false);
        ivs{jj} = vertcat(ivs{jj}{:});
    end
    
    dataTable = table(dv,ivs{:},'VariableNames',[dvNames(ii) ivNames]);
    
    model = sprintf('%s ~ %s + %s + %s + %s + (1|%s)',dvNames{ii},ivNames{:});
    
    lmes{ii} = fitlme(dataTable,model);
    glmes{ii} = fitglme(dataTable,model,'Distribution',distributions{ii});
end

%%

plotTypes = {'histogram' 'caseorder' 'fitted' 'lagged' 'probability' 'symmetry'};
models = [lmes;glmes];

for ii = 1:2
    for jj = 1:3
        figure
        
        for kk = 1:6
            subplot(2,3,kk);
            plotResiduals(models{ii,jj},plotTypes{kk});
        end
    end
end