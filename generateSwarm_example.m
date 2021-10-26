% generateSwarm_example.m

% stick all the parameters in an array
param_list = {'param1' 'param2' 'param3' 'paramN'}; % inputs to "my_function.m"

% make a command on a new line for each parameter
command_list = [];
for ii = 1:length(param_list)
    command_list = [command_list ...
        'run_my_function.sh '... % compiled "my_function.m"
        '/usr/local/matlab-compiler/v90 '... % matlab correct runtime
        param_list{ii}...
        '\n'];
end

% write the commands into a swarm file
file_handle = fopen('myjobs.swarm', 'w+'); % 'w+' will overwrite if the file already exists
fprintf(file_handle, command_list);
fclose(file_handle);



%% Example swarm file: myjobs.swarm
% run_my_function.sh /usr/local/matlab-compiler/v90 param1
% run_my_function.sh /usr/local/matlab-compiler/v90 param2
% run_my_function.sh /usr/local/matlab-compiler/v90 param3
% run_my_function.sh /usr/local/matlab-compiler/v90 paramN