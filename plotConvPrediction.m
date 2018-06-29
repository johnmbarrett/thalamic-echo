function plotConvPrediction(fig,predictions,Bfull,Afull)
% deviation=convPrediction(paras,Bfull,Afull)
% Konrad's function for convolving the upstream responses to the 
% downstream responses using a gamma distribution.
%
% Bfull: 5x5x100 upstream data array; e.g. median RSC responses to RSC stim
% Afull: 5x5x100 downstream data array; e.g. median M2 responses to RSC stim
% deviation: used to optimize the fit
% paras: e.g. [0.8965    0.0409    0.6115    1.4459    0.2]
%    paras(1): alpha (k), shape; larger-->gaussian, smaller-->exponential
%    paras(2): beta (1/theta), rate; 
%    paras(3): scale?
%    paras(4): y-offset?
%    paras(5): ?
%
% Usage: e.g.
% (1) run parametricAnalysis04
% (2) click on the fig with the 5x5 subplots of traces for the upstream
%     area of interest; e.g. RSC responses to RSC stim
% (3) run the code at the bottom of this function to generate the 'Bfull' array
% (4) do the same for the M2 responses to RSC stim to generate the 'AFull'
% (5) call this function from the command line:
% paras=[0.8726    0.0392    0.9325    1.8700    0.0840]; % AAV9 ortho
% paras=[1.0523    0.0421    0.8047    0.7867    0.0605]; % AAV1 ortho
% paras=fminsearch(@convPrediction2,paras,[],Bfull,Afull)
% 
% See also: 
% parametricAnalysis04 -- generates the figure with the 5 x 5 array 
% of subplots with blue traces, used for the fits here
% 
% 2017 may
%% Main loop ------------------------------------------------------------------
if ~(isscalar(fig) && isgraphics(fig))
    Afull = Bfull;
    Bfull = predictions;
    predictions = fig;
    figure;
else
    figure(fig);
end

x = (1:100)';
kk = 0;

for ii=1:5
    for jj=1:5
        kk = kk+1;
        
        subplot(5,5,kk);
        
        hold on;
        
        plot(x,squeeze(Bfull(ii,jj,:))/3,'Color',[0 0 0.7]);
        plot(x,squeeze(Afull(ii,jj,:)),'Color',[0.7 0 0]);
        plot(x,squeeze(predictions(ii,jj,:)),'Color',[0 0.7 0]);
    end
end