function ParSpikeSort(NeurotankPath,SavePath,chanlist,recompile_flag,blocklist,extension)

% the following gives an example of input
if nargin < 1
    NeurotankPath = '/data/NIF/projects/godlovedc/rhombus20140307/';
end
if nargin < 2
    SavePath      = '/data/NIF/projects/godlovedc/SortSpikeData3/';
end
if nargin < 3
    chanlist      = 65:128;
end
if nargin < 4
    recompile_flag = 0;
end
if nargin < 5 
    blocklist = [];
end
if nargin < 6 
    extension = 'sev';
end

warning off all
script_path = fileparts(which('PreprocessNeuroData')); % this will work from within matlab but not necessarily if the code is compiled
path(path,genpath(script_path)); % adds recursively
warning on all

% If the user has asked for it, recompile the code before running the
% swarm.  Note that another user may be tying up the license.  You can type
% licenses at the biowulf command prompt to check the status
compiler_message = 0;
if recompile_flag
    fprintf('Attempting to recompile run_PreprocessNeuroData.sh and run_remake_NeuroStruct.sh.\n')
    
    compile PreprocessNeuroData
    compile remake_NeuroStruct    
    compiler_message = 1;

end

% Run the swarm and capture the standard output to figure out what the job
% numbers end up being (for use below)
% the syntax of the swarm command in the evalc statement means...
% spawn a swarm and...
% specify the file to run is the swarm file we just created and...
% specify the queue that we want to enter (the nodes that we want to use)
% are the NIMH nodes and...
% specify that each job will need N threads (cores) so that we don't
% overload any nodes and...
% catenate a single .o (standard output) and .e (standard error) file from
% all of the jobs and...
% run the swarm with a single PBS job id to make it easier for dependencies
% below.
% but first... make the swarm file that will be needed to sort the spikes
generate_swarm_PreprocessNeuroData(NeurotankPath,SavePath,chanlist,blocklist,extension)
jobid = evalc(['!swarm -f ' script_path filesep 'swarm_PreprocessNeuroData.txt -q nimh -t 4 --singleout --jobarray']);
depend_string = sprintf('-W depend=afterany:%s',jobid); % for use below to make the next job dependent on the execution of this one
fprintf(['Spawned job number ' deblank(jobid) ' as Biowulf swarm.\n'])

% put another job into the queue that will clean up and rebuild the
% NeuroStruct after it is done being made.  Use the job numbers from the
% capture statement above to make this job contingent on the others
% executing first.
% the syntax of the swarm command in the evalc statement means...
% spawn a swarm (will just be one core here) and...
% specify the file to run is the swarm file we just created and...
% specify that we want to run our job in the NIMH queue and...
% make the job dependent on the execution of the last job so
% that it won't start until after the previous one completes.
% but first... make another file that contains the appropriate command to
% remake the NeuroStruct
generate_command_remake_NeuroStruct(NeurotankPath,SavePath,chanlist,0,script_path);
jobid = evalc(['!swarm -f ' script_path filesep 'call_run_remake_NeuroStruct.txt -q nimh ' depend_string]);
fprintf(['Submitted cleanup job number ' deblank(jobid) ' to run following swarm.\n'])

if compiler_message
   fprintf('\nIt appears that you have checked out the compiler license.\nConsider being a good citizen and logging out of this matlab\nsession if you are finished using it.\n') 
end





function compile(funname)
% funname is the name of a function on the matlab path

if strcmp(funname(end-2:end),'.m')
    funname = funname(1:end-3);
end
ismfile = 2;
if exist(funname,'file') ~=ismfile
    error(['''' funname '.m' '''' 'not found on the MATLAB search path.'])
end

try
    
    curr_dir = cd(fileparts(which(funname)));
    mcc('-m',funname)
    cd(curr_dir)    
    
catch
    
    cd(curr_dir)
    fprintf('\nWarning: Could not compile %s.m. Compiler license may be checked out.\n\n',funname);
    uinput = input(['What do you want to do?\n'...
        '[W] Wait for compiler license to become free. (May take hours.)\n'...
        '[R] Run without recompiling code (if compile was called within another function).\n'...
        '[Q] Quit.\n'],...
        's');
    bad_input = 1;
    
    while bad_input
        
        switch uinput
            
            case {'w','W'}
                
                bad_input = 0;
                
                fprintf('Waiting for compiler license to become free...\n')
                need_compiler = 1;
                
                while need_compiler
                    
                    license_info = evalc('!licenses | grep compiler');
                    free_status  = regexp(license_info,'\d');
                    free_status  = free_status(end);
                    free_status  = str2double(license_info(free_status));
                    
                    if free_status == 1
                        
                        compile(funname) % did you mean recursion?
                        need_compiler = 0;
                        
                    end
                    
                    pause(1);
                    
                end
                
            case {'r','R'}
                
                bad_input = 0;
                
            case {'q','Q'}
                
                error('Execution interrupted by user.')
                
            otherwise
                
                uinput = input(['\nPlease select an option below\n'...
                    '[W] Wait for compiler license to become free. (May take hours.)\n'...
                    '[R] Run without recompiling code (if compile was called within another function).\n'...
                    '[Q] Quit.\n'],...
                    's');
                
        end
    end    
end


