function fun = convPredictionFun(kernelFun,lower,upper,isPlot,isFitToSummedData)
% Functional version of Konrad's function for convolving the upstream 
% responses to the downstream responses using a gamma distribution,
% allowing an arbitrary kernel function to be used.
%
% kernelFun: a function handle for generating the convolution kernel, must
% take a list of parameters as its first arguments and a vector
% representing time as its second.  Alternatively, you may specify the
% string 'exp' for an exponential distribution with one parameter or
% 'gamma' for a function of the form x.^a.*e.^(-b*x) (which is not,
% strictly speaking, a gamma distribution).
%
% Returns:
% A function to be used for predicting downstream responses (e.g. in M2)
% from responses to optogenetic stimulation of an upstream area (e.g. RSC).
% The function takes the following arguments:
%
% params: a vector of parameters for the fit.  The first n parameters are
% for the kernel, the n+1th parameter is a scaling factor, the n+2th the
% minimum, and the n+3th the offset
% Bfull: 5x5x100 upstream data array; e.g. median RSC responses to RSC stim
% deviation: used to optimize the fit
% Afull: 5x5x100 downstream data array; e.g. median M2 responses to RSC stim
% 
% The function returns the deviation of the predicted response from the
% recorded response, the predicted response itself, and the kernel used to
% generate the prediction
%
% Usage: e.g.
% (1) run fitOptogeneticStimModels to generate the Afull and Bfull arrays
% (2) call this function from the command line:
% paras=[0.8726    0.0392    0.9325    1.8700    0.0840]; % AAV9 ortho
% paras=[1.0523    0.0421    0.8047    0.7867    0.0605]; % AAV1 ortho
% paras=fminsearch(convPredictionFun('exp',paras,[],Bfull,Afull) % or 'gamma' or any function handle
% 
% See also: 
% parametricAnalysis04 -- generates the figure with the 5 x 5 array 
% of subplots with blue traces, used for the fits here
% 
% 2017 may
    if ischar(kernelFun)
        switch kernelFun
            case 'exp'
                kernelFun = @(p,x) exp(-p(1)*x);
            case 'gamma'
                kernelFun = @(p,x) x.^(p(1)-1).*exp(-p(2)*x);
            otherwise
                error('Unknown kernel function %s\n',kernelFun);
        end
    end
    
    if nargin < 2 || isempty(lower)
        lower = -Inf;
    end
    
    if nargin < 3 || isempty(upper)
        upper = Inf;
    end
    
    if nargin < 4
        isPlot = false;
    elseif isPlot
        fig = figure;
    end

    if nargin < 5
        isFitToSummedData = false;
    end
    
    assert(isa(kernelFun,'function_handle'),'First argument must be a function handle'); % TODO : validate order of params and x as well?
    
    function [deviation,predictions,kernel] = convPrediction(params,Bfull,Afull)
        if any(params < lower | params > upper)
            deviation = Inf;
            predictions = NaN;
            kernel = NaN;
            return
        end
        
        x = 1:60;

        kernel = kernelFun(params,x);
        kernel = kernel/sum(kernel);

        deviation = 0;

        predictions = zeros(size(Afull));

        % TODO : introduce dim argument so arrays can be arbitrary shape
        for ii = 1:5
            for jj = 1:5
                preds = max(0,conv2(squeeze(Bfull(ii,jj,21:80)),kernel')*params(end-2)-params(end-1))+params(end);
                preds = preds(1:60);
                predictions(ii,jj,21:80) = preds;

                if isFitToSummedData
                    sumA = sum(squeeze(Afull(ii,jj,21:71)));
                    sumP = sum(preds(21:71));
                    deviation = deviation+(sumA-sumP).^2/sumA.^2;
                else
                    deviation=deviation+(sum((squeeze(Afull(ii,jj,21:80))-preds).^2))/(sum(squeeze(Afull(ii,jj,21:80).^2)));
                end
            end
        end

        if isPlot
            if ~isgraphics(fig) % restore the figure if the user closes it
                fig = figure;
            end

            plotConvPrediction(fig,predictions,Bfull,Afull);
            drawnow;
        end

        deviation = deviation/25;

        disp(deviation);
        disp(params);
    end
    
    fun = @convPrediction;
end