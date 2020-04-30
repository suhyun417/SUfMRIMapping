function [] = genSwarm_reconVoxCI(nameSubjNeural, nameSubjBOLD)

% reconVoxCI(nameSubjNeural,nameSubjBOLD, iChan)

flagRecompile = 1;

% Load the data
% dirDataHome = '/procdata/parksh'; %'/data/parks20/procdata/NeuroMRI'; % Biowulf 
% dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
% dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);
% filenameNeural = sprintf('neuraltc4computeCorrcoeffCI_%s.mat', nameSubjNeural);
% filenameBOLD = sprintf('fmritc4computeCorrcoeffCI_%s.mat', nameSubjBOLD);
% 
% load(fullfile(dirDataNeural, filenameNeural))
% load(fullfile(dirDataBOLD, filenameBOLD))

% Prep the main dimension of paramters to run parallel 
% nChan = size(neuraltc, 3);
% nVox = size(fmritc, 3);
setChan = [1:12; 13:24; 25:36; 37:48]; %1:nChan; %2:48; %1:nChan;
% setStartVox = 1:1000:81920;

script_path = '/data/parks20/analysis/NeuroMRI/_compiled/_computeCorrCI'; % where the compiled script is
compiled_script = 'run_reconVoxCI.sh';
matlab_path = '/usr/local/matlab-compiler/v91'; % matlab component runtime (MCR): version must match the version of Matlab
% v91 for R2016b, v90 for R2015b, v85 for R2015a, v84 for R2014b
mcrcachestring = 'export MCR_CACHE_ROOT=/lscratch/$SLURM_JOBID;';

if flagRecompile
    if ~isdir(script_path)
        mkdir(script_path)
    end
%     delete(fullfile(script_path, '*'))    
    mcc2('-m', '-R', '-nodisplay', '-R', '-singleCompThread', '-d', script_path, 'reconVoxCI.m')
end


% make a command on a new line for each parameter

% iTemp = size(setChan,1);
for iTemp = 1:size(setChan,1)
    command_list = [];
    for ii = 1:size(setChan,2)
        curChan = num2str(setChan(iTemp,ii));
        %     for jj = 1:length(setStartVox)
        %         curStartVox = num2str(setStartVox(jj));
        
        command_list = [command_list ...
            sprintf('%s %s %s %s %s %s\n',...
            mcrcachestring,...
            fullfile(script_path, compiled_script),...
            matlab_path,...
            nameSubjNeural, nameSubjBOLD, curChan)]; % this command will ask for a different save directory for every channel
        %     end
    end
    % end
    
    % write the commands into a swarm file
    swarmFileName = sprintf('reconVoxCI_%s%s_cell%d_%d.swarm', nameSubjNeural, nameSubjBOLD,...
        setChan(iTemp,1), setChan(iTemp,end));
    file_handle = fopen(swarmFileName, 'w+'); % fopen('computeCorrcoeffCI.swarm', 'w+'); % 'w+' will overwrite if the file already exists
    fprintf(file_handle, command_list);
    fclose(file_handle);
end
