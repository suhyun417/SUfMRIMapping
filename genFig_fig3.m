% genFig_fig3.m


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
addpath(fullfile(dirProjects, 'parksh/_toolbox/hslcolormap'))

% Set directories 
nameSubjNeural = 'Spi'; % 'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% Directory for saving figures as graphic files
dirFig = fullfile(dirProjects, 'parksh/NeuralBOLD/_labNote/_figs/');


%%
cMap = [0 0 0; 230 159 0; 86 180 233; 0 158 115; 240 228 66; 0 114 178; 213 94 0; 204 121 167]./255;
marker = {'o', '*', 'x', 's', 'd', '+', '^'};
% orderIndNewClust = [3 2 7 6 5 1 4];
% cMap_newclust = cMap(orderIndNewClust,:);

% Load data
% 1) fMRI correlation maps
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr') 
% 2) Clustering results
load(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked.mat', nameSubjNeural, nameSubjBOLD)));
% % 3) Distance matrix
% load(fullfile(dirDataNeural, sprintf('distanceMatrix_Corr%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)));
% 4) Movie-driven mask
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp');

[nx ny nz] = size(movieDrivenAmp.mask_amp1);
nVox = nx*ny*nz;

% Apply movie-driven mask to correlation matrix
moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1); % change the 3D mask to 1D
matR_SU_moviemask = matR_SU(moviemask_vec,:); % 15495 voxels


%% Fig 3a_1: Correlation matrix between voxels and cells
fig3a1=figure;
set(fig3a1, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 300 200])
imagesc(matR_SU(1:11, 1:7)') %imagesc(matR_SU(1:20, 1:10)')
axis off
rgb=hslcolormap(256, 'bc.yr', 1, [0.2 1 0.2]); colormap(rgb)
caxis([-0.18 .18])
% save
print(fig3a1, fullfile(dirFig, 'fig3a1'), '-r150', '-dtiff')

% % make blue-white-red colorbar
% cval = 0.25;
% cmin = -cval; cmax = cval;
% colornum = 256;
% colorInput = [1 0 0; 1 1 1; 0 0 1];
% oldSteps = linspace(-1, 1, length(colorInput));
% newSteps = linspace(-1, 1, colornum);
% for j=1:3 % RGB
%     newmap_all(:,j) = min(max(transpose(interp1(oldSteps, colorInput(:,j), newSteps)), 0), 1);
% end
% endPoint = round((cmax-cmin)/2/abs(cmin)*colornum);
% newmap = squeeze(newmap_all(1:endPoint, :));
% figure(gcf)
% set(gca, 'CLim', [cmin cmax])
% colormap(flipud(newmap))
% set(gca, 'TickDir', 'out')
% box off
% c=colorbar;

%% Fig 3a_2: 2-D MDS plot showing K-means clustering results
D = pdist(matR_SU_moviemask', 'euclidean');
[Y2,stress,disparities] = mdscale(D,2);

matIndClust_SU = cat(2, Clustering_moviemask.resultKMeans.SU_indCluster); %cat(2, Clustering.resultKMeans.SU_indCluster);
curK = 5; %6; %4; %5; %6; %7; %4;
[sortedClust, indSortChan]=sort(matIndClust_SU(:,curK-1));
indNewCluster = [1 2 3 4 5]; %[1 3 6 4 2 5]; %[2 4 1 3]; %[1 2 3 4 5]; %[3 4 2 5 1]; %[1 2 3 4 5 6 7]; %[4 1 6 3 5 2 7]; % cluster #4 is now cluster 1

% 
% cMap = cMap([3 4 5 2 1 6 7 8], :);
fig3a2=figure;
set(fig3a2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1300 600 700 700])
for iC=1:curK
    
    iK = indNewCluster(iC);
    figure(fig3a2); 
    curChan = indSortChan(sortedClust==iK);
    plot(Y2(curChan,1), Y2(curChan,2), 'o-','LineWidth', 2, 'MarkerSize', 10,...
        'MarkerEdgeColor','k', 'MarkerFaceColor', cMap(iC,:), 'Color', cMap(iC,:));
            text(Y2(curChan,1)+1, Y2(curChan,2), paramCorr.validChanID(curChan,:))
    hold on;
end

axis square

xlim([-25 30])
ylim([-13 13])
set(gca, 'XTick', [], 'YTick', [], 'LineWidth', 2,  'Box', 'off')
% print(fig3a2, fullfile(dirFig, 'figS4_MDS_6cluster'), '-depsc')

% save
print(fig3a2, fullfile(dirFig, 'fig3a_2'), '-r300', '-dtiff')
print(fig3a2, fullfile(dirFig, 'fig3a_2'), '-depsc')


%% Fig 3b: Explained variance elbow plot
setK = paramClustering.setK; %Clustering.setK;

matWSS=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
    indClust = Clustering_moviemask.resultKMeans(iK).SU_indCluster; %Clustering.resultKMeans(iK).SU_indCluster;
    [sortedClust, indSortedChan]=sort(indClust);
    
%     tExpVar=[];
%     for ii = 1:curK
%         tExpVar(ii,1) = Clustering.resultKMeans(iK).SU_sumD(ii)/(2*sum(sortedClust==ii));
%     end
    
    matWSS(iK,1) = sum(Clustering_moviemask.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
%     matExpVar(iK,1) = sum(tExpVar);

end

[a, c, totalSS] = kmeans(matR_SU_moviemask', 1); %kmeans(matR_SU', 1);
betweenSS = totalSS-matWSS;
% totalVar = totalSS/(2*size(matR_SU,2)); %totalD/(2*size(matR_SU,2));

propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;


fig3b=figure;
set(fig3b, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 340 340])

curK=5;
plot(setK, propExplained, 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 8); hold on
plot(curK, propExplained(curK-1), 'ko', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'MarkerSize', 8)
xlim([2 12])
ylim([0.5 0.9])
set(gca, 'TickDir', 'out', 'LineWidth', 2, 'Box', 'off', 'TickLength', [.025 .05])
set(gca, 'YTick', 0.5:.1:.9)

% save
print(fig3b, fullfile(dirFig, 'fig3b'), '-depsc')

%% 3D MDS
% [Y,stress,disparities] = mdscale(D,3);
%     
% matIndClust_SU = cat(2, Clustering_moviemask.resultKMeans.SU_indCluster); %cat(2, Clustering.resultKMeans.SU_indCluster);
% curK = 5; %6; %7; %4;
% [sortedClust, indSortChan]=sort(matIndClust_SU(:,curK-1));
% 
% fig3a2=figure;
% set(fig3a2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1300 600 700 700])
% % Using PLOT3: main MDS space
% indNewCluster = [1 2 3 4 5 6 7]; %[4 1 6 3 5 2 7]; % cluster #4 is now cluster 1
% 
% for iC=1:curK
%     
%     iK = indNewCluster(iC);
%     figure(fig3a2); 
%     curChan = indSortChan(sortedClust==iK);
%     p_org=plot3(Y(curChan,1), Y(curChan,2), Y(curChan,3), 'o-','LineWidth', 2, 'MarkerSize', 12,...
%         'Color', cMap(iC,:), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cMap(iC,:)); 
%     text(Y(curChan,1)+1, Y(curChan,2), Y(curChan,3), paramCorr.validChanID(curChan,:))
%     hold on;
% end
% 
% % Using PATCH
% % vert = Y;
% faces = padcat(indSortChan(sortedClust==1), indSortChan(sortedClust==2), indSortChan(sortedClust==3),...
%     indSortChan(sortedClust==4), indSortChan(sortedClust==5), indSortChan(sortedClust==6),...
%     indSortChan(sortedClust==7));
% faces = faces';
% % p_org = patch('Faces', faces, 'Vertices', vert, 'Marker', 'o'); hold on;
% % set(p_org, 'FaceColor', 'none', 'EdgeColor', 'flat', 'MarkerFaceColor', 'flat', 'FaceVertexCData', cdata)
% 
% % Draw projections of each dimension
% % set ideal view
% view([168.5 28])  %      view([44 44])
% % get axis limits
% axis tight
% xl = get(gca, 'XLim');
% yl = get(gca, 'YLim');
% zl = get(gca, 'ZLim');
% 
% % coordinates for background
% % x-y plane (bottom)
% bg_xy = [xl(2) xl(2) xl(1) xl(1); yl(1) yl(2) yl(2) yl(1); zl(1) zl(1) zl(1) zl(1)]; % background coords for xy plane
% vert_xy = cat(2, Y(:,1:2), repmat(zl(1), size(Y,1),1)); % data coords for xy plane
% % x-z plane (side)
% bg_xz = [xl(2) xl(2) xl(1) xl(1); yl(1) yl(1) yl(1) yl(1); zl(2) zl(1) zl(1) zl(2)]; % background coords for xy plane
% vert_xz = cat(2, Y(:,1), repmat(yl(2), size(Y,1),1), Y(:,3)); % data coords for xz plane
% % y-z plane (side)
% bg_yz = [xl(1) xl(1) xl(1) xl(1); yl(1) yl(1) yl(2) yl(2); zl(2) zl(1) zl(1) zl(2)]; % background coords for xy plane
% vert_yz = cat(2,repmat(xl(1), size(Y,1), 1), Y(:,2:3)); % data coords for yz plane
% 
% figure(fig3a2); hold on
% % grayVal = 0.95;
% % grayBG_xy = fill3(bg_xy(1,:), bg_xy(2,:), bg_xy(3,:), ones(1,3).*grayVal, 'EdgeColor', 'none'); 
% % grayBG_xz = fill3(bg_xz(1,:), bg_xz(2,:), bg_xz(3,:), ones(1,3).*grayVal, 'EdgeColor', 'none'); % for y-z plane
% % grayBG_yz = fill3(bg_yz(1,:), bg_yz(2,:), bg_yz(3,:), ones(1,3).*grayVal, 'EdgeColor', 'none'); % for y-z plane
% p(1) = patch('Faces', faces, 'Vertices', vert_xy, 'Marker', 'o'); hold on;
% % p(2) = patch('Faces', faces, 'Vertices', vert_xz, 'Marker', 'o'); hold on;
% % p(3) = patch('Faces', faces, 'Vertices', vert_yz, 'Marker', 'o'); hold on;
% 
% indOldCluster = [2 6 4 1 5 3 7]; % kmeans clustering #1 is now cell group 2 in main figure
% cdata=cMap(indOldCluster(matIndClust_SU(:,curK-1)),:);
% set(p, 'FaceColor', 'none', 'EdgeColor', 'none', 'MarkerFaceColor', 'flat', 'FaceVertexCData', cdata, 'MarkerSize', 4)
% 
% % do some extra work to make it look better
% % line for y-z plane
% l=line([bg_yz(1,:); bg_yz(1,[2 3 4 1])], [bg_yz(2,:); bg_yz(2, [2 3 4 1])], [bg_yz(3,:); bg_yz(3, [2 3 4 1])], 'Color', 'k', 'LineWidth', 2);
% set(l(1:2), 'LineStyle', ':')
% % line for x-z plane
% l=line([bg_xz(1,:); bg_xz(1,[2 3 4 1])], [bg_xz(2,:); bg_xz(2, [2 3 4 1])], [bg_xz(3,:); bg_xz(3, [2 3 4 1])], 'Color', 'k', 'LineWidth', 2);
% set(l(2:3), 'LineStyle', ':')
% % line for x-y plane
% l=line([bg_xy(1,:); bg_xy(1,[2 3 4 1])], [bg_xy(2,:); bg_xy(2, [2 3 4 1])], [bg_xy(3,:); bg_xy(3, [2 3 4 1])], 'Color', 'k', 'LineWidth', 2);
% set(l(3:4), 'LineStyle', ':')
% 
% % xlabel('axis 1'); ylabel('axis 2'); zlabel('axis 3'); % just for check
% set(gca, 'XTick', [], 'YTick', [], 'ZTick', [], 'Box', 'off')


%% 2d MDS
% [Y2,stress,disparities] = mdscale(D,2);
% 
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% plot(Y2(:,1),Y2(:,2),'o','LineWidth',2, 'MarkerSize', 8);
% text(Y2(:,1)+1, Y2(:,2), paramCorr.validChanID)
% % 
% figMDS=figure;
% set(figMDS, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 600 550 550])
% for iK = 1:8
%     
%     matIndClust_SU = cat(2, Clustering_moviemask.resultKMeans.SU_indCluster); %
%     curK = paramClustering.setK(iK); %Clustering.setK(iK);
%     [sortedClust, indSortChan]=sort(matIndClust_SU(:,curK-1));
%     
%     for ii=1:curK
%         
%         figure(figMDS);
%         curChan = indSortChan(sortedClust==ii);
%         plot(Y2(curChan,1), Y2(curChan,2), 'o-','LineWidth', 2, 'MarkerSize', 8,...
%             'MarkerEdgeColor', cMap(ii,:), 'MarkerFaceColor', cMap(ii,:), 'Color', cMap(ii,:));
% %         text(Y2(curChan,1)+1, Y2(curChan,2), paramCorr.validChanID(curChan,:))
%         hold on;
%     end
%     
%     axis square
%     title(sprintf('MDS plot: cluster # = %d', curK))
%     input('')
%     
% end

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

