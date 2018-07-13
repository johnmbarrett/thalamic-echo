topDir = 'Z:\LACIE\Manuscripts\2018 in vivo LSPS Ntsr1 etc\data';
cd(topDir);

dates = dir;
%%
dates = vertcat(dates(vertcat(dates.isdir) & arrayfun(@(s) ~isempty(regexp(s.name,'[0-9]{8}','once')),dates)).name);
nDates = size(dates,1);
%%

mps = cell(1,nDates);
xy = cell(1,nDates);
probeLocations = cell(1,nDates);
%%
close all

for ii = 1:nDates-1
    cd([topDir '\' dates(ii,:)]);
    
    mps{ii} = dir;
    
    tokens = arrayfun(@(s) regexp(s.name,'MP-10-100-([0-9]+)x([0-9]+)','tokens'),mps{ii},'UniformOutput',false);
    
    good = vertcat(mps{ii}.isdir) & ~cellfun(@isempty,tokens);
    mps{ii} = mps{ii}(good);
    tokens = tokens(good);
    
    probeLocations{ii} = zeros(numel(mps{ii}),min(4,nProbes));
    xy{ii} = zeros(numel(mps{ii}),2);
    
    for jj = 1:numel(mps{ii})
        cd([topDir '\' dates(ii,:) '\' mps{ii}(jj).name]);
        
        load([mps{ii}(jj).name '.mat']); % always right?
        
        X = str2double(tokens{jj}{1}{1});
        Y = str2double(tokens{jj}{1}{2});
        xy{ii}(jj,:) = [X Y];
        MP_ResStart_Amp_Vol_2;
        
        figs = findobj('Type','figure');
        close(figs(setdiff(1:numel(figs),(4+(4*(nProbes>4))):4:numel(figs)))); % figs array is backwards, we don't care about thalamus
        
        for kk = 1:min(4,nProbes)
            fig = figs(numel(figs)+4-4*kk);
            figure(fig);
            p = impoint;
            p = getPosition(p);
            probeLocations{ii}(jj,kk) = p(1);
            close(fig);
        end
    end
end
%%

probeLocations = cellfun(@floor,probeLocations,'UniformOutput',false);
cd(topDir);
save('probe_locations.mat','dates','mps','xy','probeLocations');