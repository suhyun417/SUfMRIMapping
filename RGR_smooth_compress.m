step4_create_fullRGR;

smoRGR = zeros(30,3360);
sqrt_RGR = zeros(30,3360);
log_RGR = zeros(30,3360);
sqrt_smoRGR = zeros(30,3360);
log_smoRGR = zeros(30,3360);

%try smooth and/or compression on the regressors
kernelLen = 8;
boxCar = ones(1,kernelLen) ./ kernelLen;

for nowNdx = 1 : 30 %
    rawData = fullRGR(nowNdx,:);
    sqrt_raw = sqrt(rawData);
    log_raw = log(rawData);
    %make sure no imagineray or inf in there
    for i = 1 : length(sqrt_raw)
        if isinf(sqrt_raw(i)) || ~isreal(sqrt_raw(i))
            sqrt_raw(i) = 0;
        end
        if isinf(log_raw(i)) || ~isreal(log_raw(i))
            log_raw(i) = 0;
        end
    end
    %clean the log/sqrt data
    
    smoData = conv( rawData, boxCar, 'same');
    smoData(1: kernelLen ./2) = rawData(1:kernelLen./2);
    smoData(end- kernelLen./2 : end) = rawData( end-kernelLen./2 : end);

    sqrt_smoData = conv(sqrt_raw, boxCar, 'same');  
    sqrt_smoData(1: kernelLen ./2) = sqrt_raw(1:kernelLen./2);
    sqrt_smoData(end- kernelLen./2 : end) = sqrt_raw( end-kernelLen./2 : end);

    log_smoData = conv(log_raw, boxCar, 'same');  
    log_smoData(1: kernelLen ./2) = log_raw(1:kernelLen./2);
    log_smoData(end- kernelLen./2 : end) = log_raw( end-kernelLen./2 : end);
    
        
%     figure;
%     hold on;
%     %plot(rawData, 'k-');
%     %plot(smoData, 'r-');
%     plot(sqrt_raw, 'b-');
%     plot(sqrt_smoData, 'g-');
%     plot(log_raw, 'c-');
%     plot(log_smoData, 'm-');
%    
%     hold off;
%     
%     zScore_raw = zscore(smoData);
%     zScore_sqrt = zscore(sqrt_smoData);
%     zScore_log = zscore(log_smoData);
%     figure;
%     hold on;
%     plot(zScore_raw, 'k-');
%     plot(zScore_log, 'm-');
%     plot(zScore_sqrt, 'b-');
%     hold off;
%     
%     
%     
%     
%     
    
    
    
smoRGR(nowNdx,:) = smoData;
sqrt_RGR(nowNdx,:) = sqrt_raw;
log_RGR(nowNdx,:) = log_raw;
sqrt_smoRGR(nowNdx,:) = sqrt_smoData;
log_smoRGR(nowNdx,:) = log_smoData;
    
    
end