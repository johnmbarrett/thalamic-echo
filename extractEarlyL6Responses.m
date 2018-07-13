topDir = 'Z:\LACIE\Manuscripts\2018 in vivo LSPS Ntsr1 etc\data';
cd(topDir);

dates = dir;
%%
dates = vertcat(dates(vertcat(dates.isdir) & arrayfun(@(s) ~isempty(regexp(s.name,'[0-9]{8}','once')),dates)).name);
nDates = size(dates,1);
%%

probeLocations = cell(1,nDates);
%%
close all

for ii = 1:nDates-1
    cd([topDir '\' dates(ii,:)]);
    
    mps = dir;
    
    tokens = arrayfun(@(s) regexp(s.name,'MP-10-100-([0-9]+)x([0-9]+)','tokens'),mps,'UniformOutput',false);
    
    good = vertcat(mps.isdir) & ~cellfun(@isempty,tokens);
    mps = mps(good);
    tokens = tokens(good);
    
    probeLocations{ii} = zeros(numel(mps),min(4,nProbes));
    
    for jj = 1:numel(mps)
        cd([topDir '\' dates(ii,:) '\' mps(jj).name]);
        
        load([mps(jj).name '.mat']); % always right?
        
        X = str2double(tokens{jj}{1}{1});
        Y = str2double(tokens{jj}{1}{2});
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
save('probe_locations.mat','dates','mps','tokens','probeLocations');