% parametric stim of bilateral S1 and thalamus
folder = 'Z:\LACIE\Manuscripts\2018 in vivo rig\data\20161010\S1-Thal\ArsBrsCrs\Par_PSTH';

[~,params,~,~,mpsth] = loadPSTH(true);

thalamicStim = params.X == 100;
% psth = psth(:,:,:,:,thalamicStim);
% cpsth = cpsth(:,:,:,thalamicStim);
% tpsth = tpsth(:,:,:,thalamicStim);
mpsth = mpsth(:,:,thalamicStim);
params = params(thalamicStim,:);

%%

thalamusPSTH = squeeze(mpsth(:,3,:));

%%

stimulusArray = zeros(size(thalamusPSTH));

for ii = 1:25
    stimulusArray(100+(1:params.PulseWidth(ii)),ii) = params.PercentPower(ii)/100;
end

%%

[~,sortIndices] = sortrows([params.PulseWidth params.PercentPower]);
Afull = reshape(thalamusPSTH(:,sortIndices)',5,5,[]);
Bfull = reshape(stimulusArray(:,sortIndices)',5,5,[]);

%%

kernelFun = convPredictionFun('exp',[],[],true);

%%

initialParams = [10 10 0 0];
[~,initialPrediction] = kernelFun(initialParams,Bfull,Afull);
%%
% figure
for nn = 1:25
    subplot(5,5,nn);
    hold on;
    [ii,jj] = ind2sub([5 5],sortIndices(nn));
    plot(squeeze(Bfull(ii,jj,:)));
    plot(squeeze(Afull(ii,jj,:)));
    plot(squeeze(initialPrediction(ii,jj,:)));
    ylim([0 1]);
end