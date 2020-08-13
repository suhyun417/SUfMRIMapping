function [] = genSwarm_clusteringMultiplePatches %(nameSubjNeural, nameSubjBOLD)

% doClusteringCorrMap_multipleSubjects_multiplePatches_parallel(stK, edK, flagSave)

flagRecompile = 1;

% Prep the main dimension of paramters to run parallel 
setStK = 2:2:40;
setEdK = setStK+1;

script_path = '/data/parks20/analysis/NeuroMRI/_compiled/_clusteringMultiplePatches'; % where the compiled script is
compiled_script = 'run_doClusteringCorrMap_multipleSubjects_multiplePatches_parallel.sh';
matlab_path = '/usr/local/matlab-compiler/v94'; % matlab component runtime (MCR): version must match the version of Matlab
mcrcachestring = 'export MCR_CACHE_ROOT=/lscratch/$SLURM_JOBID;';
flagSave = '1';

if flagRecompile
    if ~isdir(script_path)
        mkdir(script_path)
    end
%     delete(fullfile(script_path, '*'))
    
    mcc2('-m', '-R', '-nodisplay', '-R', '-singleCompThread', '-d', script_path, 'doClusteringCorrMap_multipleSubjects_multiplePatches_parallel.m')
end

% make a command on a new line for each parameter
command_list = [];
for ii = 1:length(setStK)    
    
    stK = num2str(setStK(ii));
    edK = num2str(setEdK(ii));
    
    command_list = [command_list ...
        sprintf('%s %s %s %s %s %s\n',...
        mcrcachestring,...
        fullfile(script_path, compiled_script),...
        matlab_path,...
        stK, edK, flagSave)]; % this command will ask for a different save directory for every channel
end

% write the commands into a swarm file
swarmFileName = 'clusteringMultiplePatches.swarm';
file_handle = fopen(swarmFileName, 'w+'); % fopen('computeCorrcoeffCI.swarm', 'w+'); % 'w+' will overwrite if the file already exists
fprintf(file_handle, command_list);
fclose(file_handle);
