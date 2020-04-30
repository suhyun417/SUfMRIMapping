
cd /projects/hungc/ECoG_naturalMovie_project/JulianWorkSpace/
load('HighLevelRGR.mat')
for mov = 1  : 3 
    fName = ['Movie' num2str(mov) '_10fps_rgr.mat'];
    
    rgrInfo(mov)= load(fName);
end
clear mov

lowLevelRGR_4Hz = zeros(19, 240 * 15);
for mov = 1 : 3
    for nowREG = 1 : 19
        nowRAW = rgrInfo(mov).RGR10fps.regressors(nowREG,:); %1x3000 
        nowInterp = resample( nowRAW, 2,5); %resmaple from 10hz to 4hz sampling rate
        
        lowLevelRGR_4Hz( nowREG, (1 + (mov-1).*1200) : ( mov .*1200) ) = nowInterp;
        
        
    end
end

%% Regressor storage in fullRGR variable

%1 'Luminance'    
%2 'Contrast'    
%3 'Speed'    
%4 'Rightward Motion'
%5 'Leftward Motion'    
%6 'SF (0.0 to 0.2)'    
%7 'SF (3.0 to 100.0)'
%8 'SF Ratio'    
%9 'Beta Contrast'    
%10 'Gamma Contrast'
%11 'Motion Contrast'    
%12 'Motion Dir Beta'    
%13 'Motion Dir Gamma'
%14 'Motion Div Mean'    
%15 'Motion Div Rectify'    
%16 'Motion Div STD'
%17 'Motion Div Beta'    
%18 'Motion Div Gamma'    
%19 'Scene Cuts'
%20 Faces
%21 one face
%22 body parts
%23 conspecifics
%24 humans
%25 any animal
%26 dyadic interaction
%27 One Face (full)
%28 One Face (profile)
%29 Faces (full)
%30 Faces (profile)


%construct the 56 sec x 15 clips long regressor
fullRGR = zeros(30, 224*15);
for nowREG = 1 : 30
    counter = 1;
    for clip = 1 : 15
        
            nowNdx = 9+(clip-1)* 240;
           % nowNdx
            if nowREG <= 19                
                fullRGR(nowREG,counter: counter+224-1) =  lowLevelRGR_4Hz(nowREG, nowNdx: nowNdx+224-1);
                
            elseif nowREG == 20 %face 
                nowdata = RAWRGR.rgrs(:,8) + RAWRGR.rgrs(:,9); %1x 3600entries
                fullRGR(nowREG, counter: counter+224-1) = nowdata(nowNdx: nowNdx+224-1);
            elseif nowREG == 21 %one-face
                nowdata = RAWRGR.rgrs(:,46) + RAWRGR.rgrs(:,47); %1x3600 entries
                fullRGR(nowREG, counter: counter+224-1) = nowdata(nowNdx: nowNdx+224-1);
                
            elseif nowREG == 22 %body parts
                nowdata = RAWRGR.rgrs(:,10) + RAWRGR.rgrs(:,11) + RAWRGR.rgrs(:,12);
                fullRGR(nowREG, counter: counter+224-1) = nowdata(nowNdx: nowNdx+224-1);
            elseif nowREG == 23 %conspecifics
                nowdata = RAWRGR.rgrs(:,2);
                fullRGR(nowREG, counter: counter+224-1) = nowdata(nowNdx: nowNdx+224-1);
            elseif nowREG == 24 %humans
                nowdata = RAWRGR.rgrs(:,6);
                fullRGR(nowREG, counter: counter+224-1) = nowdata(nowNdx: nowNdx+224-1);
            elseif nowREG == 25 %any animal
                nowdata = RAWRGR.rgrs(:,2) + RAWRGR.rgrs(:,3) + RAWRGR.rgrs(:,4) + RAWRGR.rgrs(:,5) + RAWRGR.rgrs(:,6);
                fullRGR(nowREG, counter: counter+224-1) = nowdata(nowNdx: nowNdx+224-1);
            elseif nowREG == 26 %dyadic interaction
                nowdata = RAWRGR.rgrs(:,33);
                fullRGR(nowREG, counter: counter+224-1) = nowdata(nowNdx: nowNdx+224-1);
            elseif nowREG == 27 %One Face (full)
                nowdata = RAWRGR.rgrs(:,46);
                fullRGR(nowREG, counter: counter+224-1) = nowdata(nowNdx: nowNdx+224-1);      
            elseif nowREG == 28 %One Face (profile)
                nowdata = RAWRGR.rgrs(:,47);
                fullRGR(nowREG, counter: counter+224-1) = nowdata(nowNdx: nowNdx+224-1);
            elseif nowREG == 29 %Faces (full)
                nowdata = RAWRGR.rgrs(:,8);
                fullRGR(nowREG, counter: counter+224-1) = nowdata(nowNdx: nowNdx+224-1); 
            elseif nowREG == 30 %Faces (profile)
                nowdata = RAWRGR.rgrs(:,9);
                fullRGR(nowREG, counter: counter+224-1) = nowdata(nowNdx: nowNdx+224-1);
            end
            counter = counter+224;
        
        
    end % for clip 
end %for nowREG

% %smooth the data
%  boxCar = ones(1,50) ./ 50;
%  s_gInfoEvenAll =  conv(gInfoEvenAll, boxCar, 'valid');
%  s_gInfoOddAll = conv(gInfoOddAll, boxCar, 'valid');
%  s_gInfoControl = conv(gInfoControl, boxCar, 'valid');
%  for reg = 1 : 22
%      s_fullRGR(reg,:) = conv(fullRGR(reg,:), boxCar, 'valid');
%  end
% clear corrMap;
% for reg = 1 : 22
%    corrMap(reg,1) = corr(s_gInfoEvenAll, s_fullRGR(reg,:)');
%    corrMap(reg,2) = corr(s_gInfoOddAll, s_fullRGR(reg,:)');
%    corrMap(reg,3) = corr(s_gInfoControl, s_fullRGR(reg,:)');
% end
% figure; imagesc(corrMap);
% caxis([-0.4 0.4]);
% 
% figure; 
% hold on;
% plot(s_gInfoEvenAll);
% plot(s_gInfoOddAll, 'r-');
% plot(s_fullRGR(9,:), 'k-');
% hold off;



