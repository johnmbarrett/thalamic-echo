topDir = 'Z:\LACIE\Manuscripts\2018 in vivo LSPS Ntsr1 etc\data';
cd(topDir);

load('probe_locations.mat');
%%

masterFile = 'all_experiments.xlsx';

if exist(masterFile,'file')
%     error('Master spreadsheet already exists. Go in and comment out this line of code if you really want to do this.');
end

%%

nRecordings = sum(cellfun(@numel,mps(1:numel(probeLocations))));

masterTable = table(repmat(datetime('now','TimeZone','America/Chicago'),nRecordings,1),cell(nRecordings,1),zeros(nRecordings,1),zeros(nRecordings,1),cell(nRecordings,1),'VariableNames',{'Date' 'MPFolder' 'X' 'Y' 'ProbeLocations'});

%%

nn = 0;

for ii = 1:numel(probeLocations)
    for jj = 1:numel(mps{ii})
        nn = nn + 1;
        
        masterTable.Date(nn) = datetime(dates(ii,:),'InputFormat','yyyyMMdd','TimeZone','America/Chicago');
        masterTable.MPFolder{nn} = mps{ii}(jj).name;
        masterTable.X(nn) = xy{ii}(jj,1);
        masterTable.Y(nn) = xy{ii}(jj,2);
        masterTable.ProbeLocations{nn} = probeLocations{ii}(jj,:);
    end
end

%%

writetable(masterTable,masterFile);