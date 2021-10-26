function generate_command_remake_NeuroStruct(NeurotankPath,SavePath,chanlist,MRI_flag,script_path)
% This function generates and saves the command that will be needed to call
% run_remake_NeuroStruct.sh will all appropriate inputs.  The command can
% then be called by the pbs batch system on biowulf

if nargin < 5, script_path = '/data/NIF/projects/godlovedc/Code/SortSpikes'; end
if nargin < 4, MRI_flag = []; end

compiled_script = 'run_remake_NeuroStruct.sh'; % the name of the compiled mat code that will be called in the command
bin_path        = '/usr/local/matlab64'; % location of matlab
chanstr         = ['"' mat2str(chanlist) '"'];
MRI_flag        = num2str(MRI_flag);

command_list = sprintf('%s%s%s %s %s %s %s %s',script_path,filesep,compiled_script,bin_path,NeurotankPath,SavePath,chanstr,MRI_flag); % this command will ask for a different save directory for every channel

h = fopen([script_path filesep 'call_run_remake_NeuroStruct.txt'],'w+');
fprintf(h,command_list);
fclose(h);
