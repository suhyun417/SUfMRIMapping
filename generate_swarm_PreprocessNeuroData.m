function generate_swarm_PreprocessNeuroData(NeurotankPath,SavePath,chanlist,blocklist,extension) 
% This function makes the swarm file that will be needed to spawn a swarm 
% that will sort each channel on a different core.  A swarm file is just a 
% list of commands with a single command per line.  Each command is spawned 
% as an individual process in the swarm. 

script_path     = [fileparts(which('PreprocessNeuroData')) filesep]; % where the compiled mat code is located and where the swarm file will be saved 
compiled_script = 'run_PreprocessNeuroData.sh'; % the name of the compiled mat code that will be called by the swarm 
bin_path        = '/usr/local/matlab64'; % location of matlab (on every node in the swarm) 
block_str       = []; 
for ii = 1:length(blocklist) 
    block_str   = [block_str '''' char(blocklist(ii)) '''' ';']; 
end 
block_str       = ['"{' block_str(1:end-1) '}"']; % this foolishness is b/c 
                                                  % I only know how to pass 
                                                  % strings as arguments in 
                                                  % bash :-( 


command_list = []; 
for ii = 1:length(chanlist) % this assumes that channel is the best dimension for a swarm 
                            % if muliple files are being run, it might be 
                            % good to run the swarm across that dimension 
                            % as well. 
     
    curr_channel = num2str(chanlist(ii));     
     
    curr_command = sprintf('%s%s %s %s %s%s/ %s %s %s %s',script_path,... 
        compiled_script,... 
        bin_path,... 
        NeurotankPath,... 
        SavePath,... 
        curr_channel,... 
        curr_channel,... 
        block_str,... 
        extension,... 
        script_path); % this command will ask for a different save directory for every channel     
     
    command_list = [command_list curr_command '\n']; 
     
end 

h = fopen([script_path 'swarm_PreprocessNeuroData.txt'],'w+'); 
fprintf(h,command_list); 
fclose(h); 
