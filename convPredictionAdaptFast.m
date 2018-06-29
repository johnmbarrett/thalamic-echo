function [deviation,predictions,pastFilter]=convPredictionAdaptFast(paras,Bfull,Afull)
% deviation=convPrediction(paras,Bfull,Afull)
% Konrad's second function, this one for convolving the photostimulus pulse itself to the 
% local responses.
%
% Bfull: 5x5x100 stimulus array -- 
% NOTE: potential source of confusion: you can generate this as the "C" array using the
% code at the bottom of this function; but note that within this function
% it C is then assigned to Bfull, and Bfull is assigned to Afull (for
% historical reasons ...)
% Afull: 5x5x100 upstream data array; e.g. median RSC responses to RSC stim
% deviation: used to optimize the fit
% paras: e.g. 1.0e+03 *[0.1517    1.3090    0.0000   -0.0035    0.0007]
%
% Usage: e.g.
% (1) run parametricAnalysis04
% (2) click on the fig with the 5x5 subplots of traces for the upstream
%     area of interest; e.g. RSC responses to RSC stim
% (3) run the code at the bottom of this function to generate the 'C' array
% in the workspace
% (4) do the same for the M2 responses to RSC stim to generate the 'BFull'
% in the workspace
% (5) call this function from the command line:
% paras=1.0e+03*[1.75    0.0   -0.0035    0.0007    0.025]; % AAV9 ortho
% paras=1.0e+03*[ 3000    0.0   -0.0035    0.00007    0.07]; % AAV1 ortho
% paras=fminsearch(@convPredictionAdapt,paras,[],C, Bfull)
% 
% See also: 
% parametricAnalysis04 -- generates the figure with the 5 x 5 array 
%   of subplots with blue traces, used for the fits here
% convPrediction2 -- similar to this, but for fitting the downstream
% responses
% 
% 2017 may
% ------------------------------------------------------------------
deviation=0;
predictions = zeros(size(Bfull));
for i=1:5
    for j=1:5
        preds=squeeze(Bfull(i,j,:));
        pastFilter=ones(60,1);
        pastFilter=pastFilter/sum(pastFilter);
        normalizer=conv(preds,pastFilter);
        normalizer=normalizer(1:length(preds));
        preds=preds./abs(normalizer+paras(2));
        preds=interp1(1:length(preds),preds,(1:length(preds))+paras(3))';
        preds(isnan(preds))=0;
        preds=preds*paras(4);
        predictions(i,j,:) = preds;
        deviation=deviation+sqrt(sum((squeeze(Afull(i,j,:))-preds).^2));
    end
end
disp(deviation)
disp(paras)