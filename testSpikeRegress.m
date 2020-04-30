
clear;

load TorRhoSigSpi_movieTS_SU_movie123_30fps.mat

testChan = [matTS{1}(1:9000,:); matTS{2}(1:9000,:); matTS{3}(1:9000,:)];

m                 = matfile('combineMovie_PCA.mat');
PCscoreForRegress = m.PCscoreForRegress(:,:);
pcaPredict        = zeros(size(testChan,1)/3, size(testChan,2));
reg_st            = zeros(size(testChan,2), 4);
reg_coef          = zeros(1+size(PCscoreForRegress,2),size(testChan,2));

% m_diff                 = matfile('combineMovie_diff_PCA.mat');
% PCscoreForRegress_diff = m_diff.PCscoreForRegress_diff;
% pcaDiff_Predict        = zeros(size(testChan,1)/3, size(testChan,2));
% reg_diff_st            = zeros(size(testChan,2), 4);
% reg_diff_coef          = zeros(1+size(PCscoreForRegress_diff,2),size(testChan,2));


for ChList = 1:135
    disp([mfilename ' >> processing neuron #' int2str(ChList)]);
    currCh     = testChan(:,ChList);
    currCh     = reshape(currCh, 3, []); %%mean(reshape(currCh, 3,[]))';
    currCh     = sum(currCh)';
    
    [reg_coef(:,ChList),~, ~, ~, reg_st(ChList,:)] = regress(currCh, [ones(9000,1) PCscoreForRegress]);
    pcaPredict(:, ChList) = [ones(9000,1) PCscoreForRegress]*reg_coef(:, ChList);
%     st
%     plot(1:9000,currCh, 'r.-', 1:9000,pcaPredict(:, ChList),'bs-'); 
%     title(['this is neuron #' int2str(ChList) ' combined moive R^2=' num2str(reg_st(ChList,))]);

%     randPC   = PCscoreForRegress(:);
%     randPC   = randPC(randperm(length(randPC)));
%     randPC   = reshape(randPC, size(PCscoreForRegress,1),[]);
%     [reg_diff_coef(:,ChList),~, ~, ~, reg_diff_st(ChList,:)] = regress(currCh, [ones(9000,1) PCscoreForRegress_diff]);
%     pcaDiff_Predict(:, ChList) = [ones(9000,1) PCscoreForRegress_diff]*reg_diff_coef(:,ChList);
%     
% %     rand_st
%     subplot(212), plot(1:9000,currCh, 'r.-', 1:9000,pcaDiff_Predict(:,ChList),'bs-');
%     title([' combined movie with differnece R^2=' num2str(reg_diff_st(ChList,:))]);
%     drawnow; 
end    

tmp1 = corrcoef(testChan);
tmp2 = corrcoef(reg_coef(2:end,:));
corrcoef(tmp1(:),tmp2(:))

subplot(411), imagesc(tmp1);
subplot(412), imagesc(tmp2);
subplot(413), hist(reg_st(:,1),100);
subplot(414), hist(-log10(reg_st(:,3)),100);

tmpP = -log10(reg_st(:,3));


