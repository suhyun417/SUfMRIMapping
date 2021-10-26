% createMovieRGR_size_indMov.m

% DM's scene-based features (face, torso, arms, legs)
for iMovie = 1:3
    load(sprintf('/procdata/parksh/MovieRegressors/annotationMovie%d.mat', iMovie))
    sceneInfo = cat(2, sta', sto');
    nScene = length(epoch);
    validFrame_range = [1 9000]; % [(5*(iMovie-1)*60*30)+1, (5*iMovie*60*30)];
  
    matSizeRGR = [];
    for iS = 1:nScene
        matSizeRGR(sta(iS):sto(iS),1) = epoch(iS).notes.face.A;
        matSizeRGR(sta(iS):sto(iS),2) = epoch(iS).notes.torso.A;
        matSizeRGR(sta(iS):sto(iS),3) = epoch(iS).notes.arms.A;
        matSizeRGR(sta(iS):sto(iS),4) = epoch(iS).notes.legs.A;
        matSizeRGR(sta(iS):sto(iS),5) = epoch(iS).notes.viewAngle;
    end
    
    if sto(nScene) < validFrame_range(2) % add NaN to the end of the regressor
        matPad = NaN(validFrame_range(2) - sto(nScene), 5);
        matSizeRGR = cat(1, matSizeRGR, matPad);
    end
    
    sizeRGR(iMovie).face = matSizeRGR(validFrame_range(1):validFrame_range(2),1);
    sizeRGR(iMovie).torso = matSizeRGR(validFrame_range(1):validFrame_range(2),2);
    sizeRGR(iMovie).arms = matSizeRGR(validFrame_range(1):validFrame_range(2),3);
    sizeRGR(iMovie).legs = matSizeRGR(validFrame_range(1):validFrame_range(2),4);
    sizeRGR(iMovie).viewAngle = matSizeRGR(validFrame_range(1):validFrame_range(2),5);
        
    
end
    


    
%     validFrame_range = [(5*(iMovie-1)*60*30)+1, (5*iMovie*60*30)];
%     validLocHigh = [];
%     validLocHigh = locHigh(locHigh>validFrame_range(1) & locHigh<validFrame_range(2));
%     
%     % should be done frame-by-frame
%     tempS = [];
%     for iFrame = 1:length(validLocHigh)
%         indFrame = validLocHigh(iFrame)-(5*(iMovie-1)*60*30);
%         indScene = find(sum([sta<=indFrame; sto>=indFrame], 1)>1);
%         
%         if isempty(indScene)
%             continue;
%         end
%         
%         tempS(iFrame, 1) = epoch(indScene).notes.face.A;
%         tempS(iFrame, 2) = epoch(indScene).notes.torso.A;
%         tempS(iFrame, 3)= epoch(indScene).notes.arms.A;
%         tempS(iFrame, 4) = epoch(indScene).notes.legs.A;
%         %             tCheck(iFrame, 1) = indFrame;
%         %             tCheck(iFrame, 2) = indScene;
%     end

