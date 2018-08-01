nProbes = size(MP,2)/32;
Probes = 1:nProbes;

Response_start_2 = zeros(X,Y,nProbes);
Response_peak_2 = zeros(X,Y,nProbes);
Response_vol_2 = zeros(X,Y,nProbes);

for h = 1:nProbes
    Nchannels = 32;
    PreSti = 20;
    PostSti = 100;
%     X = 12;
%     Y = 12;
    Probe = Probes(h); %1:5;           %%%%%%%%%%
    channels = 1:32;     %%%%%%%%%%%
    Time_range = 1:120;
    Channel_range = nan(0,0);
    for i = 1:size(Probe,2)
    Channel_range = cat(2,Channel_range,(Probe(i)-1)*Nchannels+channels);
    end
    PSTH_p = MP(Time_range,Channel_range,:);

    %PSTH_t(:,:) = mean(PSTH,2);
    PSTH_t1 = mean(PSTH_p,2);
    PSTH_t = PSTH_t1(:,:);
%     figure; this is pre artifact subtraction, we don't care about that
%     imagesc(PSTH_t1);
    %%
%     PSTH_t(21:38,:) = repmat(mean(PSTH_t(1:19,:)),18,1);
    PSTH_t(21,:) = mean(PSTH_t(1:19,:));
    PSTH_t(22,:) = mean(PSTH_t(1:19,:));
    PSTH_t(23,:) = mean(PSTH_t(1:19,:));
    PSTH_t(24,:) = mean(PSTH_t(1:19,:));

    PSTH_t(31,:) = (PSTH_t(30,:)+PSTH_t(34,:))/2;
    PSTH_t(32,:) = (PSTH_t(31,:)+PSTH_t(34,:))/2;
    PSTH_t(33,:) = (PSTH_t(32,:)+PSTH_t(34,:))/2;

    figure;
    imagesc(PSTH_t);
    % plot(PSTH_t(:,:));
    % %ylim([0 max_y1]);
    % xlim([1 size(PSTH_t,1)]);
    % set(gca,'xtick',1:20:500,'xticklabel',-PreSti:20:500-PreSti);

    Vol_all = sum(PSTH_t(21:120,:),1);

    PSTH_t_all = sum(PSTH_t,2);
    PSTH_t_alls = smooth(PSTH_t_all,5);
    Resp_t_thd = mean(PSTH_t_alls(1:PreSti)) + 0.1 * (max(PSTH_t_alls) - mean(PSTH_t_alls(1:PreSti)));  %%%%%%%%%%%%%%%

    % JB 2018-08-01 - if the response is week, this would find the last
    % window in which the baseline creeps above the threshold, not the
    % actual response
%     for i = 1:size(PSTH_t_alls,1)-1
%         if PSTH_t_alls(i, 1) <= Resp_t_thd && PSTH_t_alls(i+1, 1) > Resp_t_thd
%             Resp_t(1,1) = i;
%         end
%         if PSTH_t_alls(i, 1) > Resp_t_thd && PSTH_t_alls(i+1, 1) <= Resp_t_thd
%             Resp_t(2,1) = i;
%         end
%     end

    [~,maxTime] = max(PSTH_t_alls((PreSti+1):end));
    maxTime = maxTime + PreSti;
    Resp_t = [0;0];
    Resp_t(1) = max(1,sum(find(PSTH_t_alls((PreSti+1):maxTime) >= Resp_t_thd,1)))+PreSti; % first above threshold sample, sum of scalar X = X, sum of [] = 0, hence wrapped in sum
    Resp_t(2) = min(PostSti,sum(find(PSTH_t_alls((maxTime+1):end) < Resp_t_thd,1))+maxTime-1);

    Resp_R = zeros(size(PSTH_t,2),3);
    PSTH_tns = zeros(size(PSTH_t,1),size(PSTH_t,2));
    for i = 1:size(PSTH_t,2)    
        Resp_R(i,1) = mean(PSTH_t(1:PreSti,i));
        Resp_R(i,2) = mean(PSTH_t(Resp_t(1,1):Resp_t(2,1),i));
        Resp_R(i,3) = (Resp_R(i,2) - Resp_R(i,1)) / Resp_R(i,1);
        if Resp_R(i,3) > 2           %%%%%%%%%%%%%%%%%%%%%%%%%%
            PSTH_tn = PSTH_t(:,i) / max(PSTH_t(Resp_t(1,1):Resp_t(2,1),i));
        else
            PSTH_tn = zeros(size(PSTH_t,1),1);
        end
        PSTH_tns(:,i) = smooth(PSTH_tn,5);
    end

    figure;
    plot(PSTH_tns(:,:));
    ylim([0 1]);
    xlim([1 size(PSTH_t,1)]);
    set(gca,'xtick',1:50:500,'xticklabel',-PreSti:50:500-PreSti);
    %axis tight


    figure;
    plot(PSTH_t_all);
    %plot(PSTH_tns(:,n));
    %ylim([0 1]);
    xlim([1 size(PSTH_t,1)]);
    set(gca,'xtick',1:50:500,'xticklabel',-PreSti:50:500-PreSti);
    %axis tight

    Resp_s = zeros(size(PSTH_tns,2),1);
    for j = 1:size(PSTH_tns,2)
        Resp_s_thd = 0.5 * max(PSTH_tns(:,j));    %%%%%%%%%%%%%%%%
        if PSTH_tns(:,j) == 0
            Resp_s(j,1) = 0;
        else
            ct = 0;
        for i = 1:size(PSTH_tns(:,j),1)-1
            if PSTH_tns(i,j) <= Resp_s_thd && PSTH_tns(i+1,j) > Resp_s_thd
                if ct == 0
                Resp_s(j,1) = i;
                ct = 1;
                end
            end
            %     if PSTH_tns(i,j) > Resp_s_thd && PSTH_tns(i+1,j) <= Resp_s_thd
            %         Resp_t(j,1) = i;
            %     end
        end
        end
    end
    Response_start = Resp_s - PreSti;
    Response_start_1 = Response_start;
    for i = 1:size(Response_start,1)
        if Response_start(i) < 0
            Response_start_1(i) = max(Response_start)+5;
        else
            Response_start_1(i) = Response_start(i);
        end

    end

    Response_start_2(:,:,h) = reshape(Response_start_1,Y,X)';

    figure;
    subplot(1,3,1);
    imagesc(Response_start_2(:,:,h));
    %contour(Response_start_2');
    axis square
    colorbar
    colormap hot
    title('Response Start Time');
    %bar(Response_start);

    Response_peak = zeros(size(PSTH_t,2),1);
    for i = 1:size(PSTH_t,2)       
        PSTH_tr = PSTH_t(:,i);% - Resp_R(i,1)) / Resp_R(i,1);   
        PSTH_trs = smooth(PSTH_tr,5);
        Response_peak(i,1) = sum(max(PSTH_trs(Resp_t(1,1):Resp_t(2,1))));
    end

    Response_peak_2(:,:,h) = reshape(Response_peak,Y,X)';
    subplot(1,3,2);
    imagesc(Response_peak_2(:,:,h));
    %contour(Response_peak_2');
    axis square
    colorbar
    colormap hot
    title('Response Peak');
    %bar(Response_peak);

    Response_vol_2(:,:,h) = reshape(Vol_all,Y,X)';
    subplot(1,3,3);
    imagesc(Response_vol_2(:,:,h));
    %contour(Response_peak_2');
    axis square
    colorbar
    colormap hot
    title('Response Volume');
%bar(Response_peak);
end

save('response_params.mat','Response_start_2','Response_peak_2','Response_vol_2');