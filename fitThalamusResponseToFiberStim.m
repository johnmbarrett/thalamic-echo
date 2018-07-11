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
Afull = reshape(thalamusPSTH(2:end,sortIndices)',5,5,[]);
Bfull = reshape(stimulusArray(2:end,sortIndices)',5,5,[]);

isPlot = false;
kernelFun = @convPredictionAdaptFast; %convPredictionFun('exp',0,Inf,isPlot,false,100,100);
%%
initialParams = [0 -10 0.001];
[~,initialPrediction] = kernelFun(initialParams,Bfull,Afull);
plotConvPrediction(initialPrediction,Bfull,Afull);
%%
optimisedParams = fminsearch(@(p) kernelFun(p,Bfull,Afull),initialParams);
[~,optimisedPrediction] = kernelFun(optimisedParams,Bfull,Afull);
plotConvPrediction(optimisedPrediction,Bfull,Afull);