% doMDS.m
%
% perform multidimensional scaling on map data

%% Settings
ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/PROJECTS';
    dirProcdata = '/Volumes/PROCDATA';
    dirLibrary = '/Volumes/LIBRARY';
else % on virtual machine
    dirProjects = '/projects';
    dirProcdata = '/procdata';
    dirLibrary = '/library';
end
    
% Add necessary toolbox % Should be 
addpath(fullfile(dirLibrary, 'matlab_utils')) % for convolution
addpath(fullfile(dirProjects, 'parksh/_toolbox/'))

% Set directories 
nameSubjNeural = 'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% addpath('/library/matlab_utils/')

cMap = [0 0 0; 230 159 0; 86 180 233; 0 158 115; 240 228 66; 0 114 178; 213 94 0; 204 121 167]./255;
marker = {'o', '*', 'x', 's', 'd', '+', '^'};
% cMap = jet(7);

% % Load the data
% nameSubjNeural = 'Tor';
% nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
% dirDataHome = '/procdata/parksh/';
% dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
% dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% % 1) fMRI correlation maps
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr') 
% 2) Clustering results
load(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)));

matIndClust_SU = cat(2, Clustering.resultKMeans.SU_indCluster);


% D = pdist(matR_SU'); 
% matD = squareform(D);

% Load distance matrix
load(fullfile(dirDataNeural, sprintf('distanceMatrix_Corr%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)));


% [Y2,stress,disparities] = mdscale(D,2);
% 
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% plot(Y2(:,1),Y2(:,2),'o','LineWidth',2, 'MarkerSize', 8);
% text(Y2(:,1)+1, Y2(:,2), paramCorr.validChanID)
% 
% figMDS=figure;
% set(figMDS, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1300 600 550 550])
% for iK = 1:8
%     
%     curK = Clustering.setK(iK);
%     [sortedClust, indSortChan]=sort(matIndClust_SU(:,curK-1));
%     
%     for ii=1:curK
%         
%         figure(figMDS);
%         curChan = indSortChan(sortedClust==ii);
%         plot(Y2(curChan,1), Y2(curChan,2), 'o','LineWidth', 2, 'MarkerSize', 8, 'MarkerEdgeColor', cMap(ii,:), 'MarkerFaceColor', cMap(ii,:));
%         text(Y2(curChan,1)+1, Y2(curChan,2), paramCorr.validChanID(curChan,:))
%         hold on;
%     end
%     
%     axis square
%     title(sprintf('MDS plot: cluster # = %d', curK))
%     input('')
%     
% end

%% 3-D MDS plot
[Y,stress,disparities] = mdscale(D,3);
figMDS3d=figure;
set(figMDS3d, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1300 600 700 700])
% for iK = 1:8
    
curK = 7; %Clustering.setK(iK);
[sortedClust, indSortChan]=sort(matIndClust_SU(:,curK-1));

cdata=cMap(matIndClust_SU(:,curK-1),:);

% Using PLOT3: main MDS space
for ii=1:curK
    
    figure(figMDS3d);
    curChan = indSortChan(sortedClust==ii);
    p_org=plot3(Y(curChan,1), Y(curChan,2), Y(curChan,3), 'o-','LineWidth', 2, 'MarkerSize', 10,...
        'Color', cMap(ii,:), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cMap(ii,:));
%     text(Y(curChan,1)+1, Y(curChan,2), Y(curChan,3), paramCorr.validChanID(curChan,:))
    hold on;
end
% Using PATCH
% vert = Y;
% faces = padcat(indSortChan(sortedClust==1), indSortChan(sortedClust==2), indSortChan(sortedClust==3),...
%     indSortChan(sortedClust==4), indSortChan(sortedClust==5), indSortChan(sortedClust==6),...
%     indSortChan(sortedClust==7));
% faces = faces';
% p_org = patch('Faces', faces, 'Vertices', vert, 'Marker', 'o'); hold on;
% set(p_org, 'FaceColor', 'none', 'EdgeColor', 'flat', 'MarkerFaceColor', 'flat', 'FaceVertexCData', cdata)

% Draw projections of each dimension
view([44 44])
% get axis limits
xl = get(gca, 'XLim');
yl = get(gca, 'YLim');
zl = get(gca, 'ZLim');

% coordinates for background
% x-y plane (bottom)
bg_xy = [xl(1) xl(2) xl(2) xl(1); yl(1) yl(1) yl(2) yl(2); zl(1) zl(1) zl(1) zl(1)]; % background coords for xy plane
vert_xy = cat(2, Y(:,1:2), repmat(zl(1), size(Y,1),1)); % data coords for xy plane
% x-z plane (side)
bg_xz = [xl(1) xl(1) xl(2) xl(2); yl(2) yl(2) yl(2) yl(2); zl(2) zl(1) zl(1) zl(2)]; % background coords for xy plane
vert_xz = cat(2, Y(:,1), repmat(yl(2), size(Y,1),1), Y(:,3)); % data coords for xz plane
% y-z plane (side)
bg_yz = [xl(1) xl(1) xl(1) xl(1); yl(1) yl(1) yl(2) yl(2); zl(2) zl(1) zl(1) zl(2)]; % background coords for xy plane
vert_yz = cat(2,repmat(xl(1), size(Y,1), 1), Y(:,2:3)); % data coords for yz plane

figure(figMDS3d); hold on
grayVal = 0.95;
grayBG_xy = fill3(bg_xy(1,:), bg_xy(2,:), bg_xy(3,:), ones(1,3).*grayVal, 'EdgeColor', 'none'); 
grayBG_xz = fill3(bg_xz(1,:), bg_xz(2,:), bg_xz(3,:), ones(1,3).*grayVal, 'EdgeColor', 'none'); % for y-z plane
grayBG_yz = fill3(bg_yz(1,:), bg_yz(2,:), bg_yz(3,:), ones(1,3).*grayVal, 'EdgeColor', 'none'); % for y-z plane
p(1) = patch('Faces', faces, 'Vertices', vert_xy, 'Marker', 'o'); hold on;
p(2) = patch('Faces', faces, 'Vertices', vert_xz, 'Marker', 'o'); hold on;
p(3) = patch('Faces', faces, 'Vertices', vert_yz, 'Marker', 'o'); hold on;
set(p, 'FaceColor', 'none', 'EdgeColor', 'none', 'MarkerFaceColor', 'flat', 'FaceVertexCData', cdata)

% do some extra work to make it look better
% line for y-z plane
l=line([bg_yz(1,:); bg_yz(1,[2 3 4 1])], [bg_yz(2,:); bg_yz(2, [2 3 4 1])], [bg_yz(3,:); bg_yz(3, [2 3 4 1])]);
[l(2:3).LineStyle] = deal(':');
[l.Color] = deal([0 0 0]);
% line for x-z plane
l=line([bg_xz(1,:); bg_xz(1,[2 3 4 1])], [bg_xz(2,:); bg_xz(2, [2 3 4 1])], [bg_xz(3,:); bg_xz(3, [2 3 4 1])]);
[l(1:2).LineStyle] = deal(':');
[l.Color] = deal([0 0 0]);
% line for x-y plane
l=line([bg_xy(1,:); bg_xy(1,[2 3 4 1])], [bg_xy(2,:); bg_xy(2, [2 3 4 1])], [bg_xy(3,:); bg_xy(3, [2 3 4 1])]);
[l(3:4).LineStyle] = deal(':');
[l.Color] = deal([0 0 0]);

xlabel('axis 1'); ylabel('axis 2'); zlabel('axis 3'); % just for check
axis off

%% save figure





% alpha = 0.7;
% pm = cat(1, p.MarkerHandle); % get hidden marker handle to control transparency of vertex markers
% faceCData = uint8(cMap(sortedClust,:).*255)';
% faceCData(4,:) = uint8(ones(1, size(faceCData,2)).*alpha*255);
% 
% [pm.FaceColorData] = deal(faceCData);
% [pm.Size] = deal(8);



% for ii=1:curK
%     
%     figure(figMDS3d);
%     curChan = indSortChan(sortedClust==ii);
%     p_org=plot3(Y(curChan,1), Y(curChan,2), Y(curChan,3), 'o-','LineWidth', 2, 'MarkerSize', 10,...
%         'Color', cMap(ii,:), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cMap(ii,:));
% %     text(Y(curChan,1)+1, Y(curChan,2), Y(curChan,3), paramCorr.validChanID(curChan,:))
%     hold on;
% end
% 
% % set a perfect view before you get the limit of each axis
% view([40 30])
% 
% % get limits
% xl = get(gca, 'XLim');
% yl = get(gca, 'YLim');
% zl = get(gca, 'ZLim');
% 
% yy = [yl(1) yl(1) yl(2) yl(2)];
% xx = [xl(2) xl(2) xl(2) xl(2)];
% zz = [zl(2) zl(1) zl(1) zl(2)];
% ff = fill3(xx, yy, zz, ones(1,3).*0.95, 'EdgeColor', 'none'); % for y-z plane
% hold on;


% %% using plot
% l = size(Y, 1);
%     
% % draw other projections
% alpha=0.5;
% set(gcf, 'Renderer', 'OpenGL')
% % x-y plane
% for ii=1:curK    
%     figure(figMDS3d);
%     curChan = indSortChan(sortedClust==ii);
%     p_xy=plot3(Y(curChan,1), Y(curChan,2), repmat(zl(1), length(curChan), 1), 'o','LineWidth', 2, 'MarkerSize', 8,...
%         'Color', cMap(ii,:), 'MarkerEdgeColor', cMap(ii,:), 'MarkerFaceColor', cMap(ii,:));
% %     p_xy.MarkerHandle.FaceColorData(4,:) = uint8(ones(1, size(p_xy.MarkerHandle.FaceColorData,2)).*alpha*255);
%     hold on;
% end
% % x-z plane
% for ii=1:curK    
%     figure(figMDS3d);
%     curChan = indSortChan(sortedClust==ii);
%     plot3(Y(curChan,1), repmat(yl(2), length(curChan), 1), Y(curChan,3), 'o-','LineWidth', 2, 'MarkerSize', 8,...
%         'Color', cMap(ii,:), 'MarkerEdgeColor', cMap(ii,:), 'MarkerFaceColor', cMap(ii,:));
%     hold on;
% end
% % y-z plane
% for ii=1:curK    
%     figure(figMDS3d);
%     curChan = indSortChan(sortedClust==ii);
%     p=plot3(repmat(xl(2), length(curChan), 1), Y(curChan,2), Y(curChan,3), 'o','LineWidth', 2, 'MarkerSize', 8,...
%         'Color', cMap(ii,:), 'MarkerEdgeColor', cMap(ii,:), 'MarkerFaceColor', cMap(ii,:));
%     hold on;
% end
%     
%     axis square
%     grid on
%     title(sprintf('MDS plot: cluster # = %d', curK))
% 





% figure(figMDS3d); hold on
% grayVal = 0.95;
% grayBG_xy = fill3(bg_xy(1,:), bg_xy(2,:), bg_xy(3,:), ones(1,3).*grayVal, 'EdgeColor', 'none'); 
% % x-y plane
% for ii=1:curK    
%     figure(figMDS3d);
%     curChan = indSortChan(sortedClust==ii);
%     p_xy=plot3(Y(curChan,1), Y(curChan,2), repmat(zl(1), length(curChan), 1), 'o','LineWidth', 2, 'MarkerSize', 8,...
%         'Color', cMap(ii,:), 'MarkerEdgeColor', cMap(ii,:), 'MarkerFaceColor', cMap(ii,:));
% %     p_xy.MarkerHandle.FaceColorData(4,:) = uint8(ones(1, size(p_xy.MarkerHandle.FaceColorData,2)).*alpha*255);
%     hold on;
% end




%     input('')
%     
% end

% for iK=1:7
% curChanInd = indSortChan(sortedClust==iK);
% plot3(Y(curChanInd,1),Y(curChanInd,2), Y(curChanInd,3),'o','LineWidth',2, 'MarkerSize', 8, ...
%     'MarkerEdgeColor', cMap(iK,:), 'MarkerFaceColor', cMap(iK,:)); %, ...
% %     'Marker', marker{iK});
% hold on
% % text(Y(curChanInd,1)+1, Y(curChanInd,2), Y(curChanInd,3), paramCorr.validChanID(curChanInd,:))
% end

