function NeuroStruct = remake_NeuroStruct(NeurotankPath,SavePath,chanlist,MRI_flag)
% The way that PreporcessNeuroData.m is called by the swarm spreads the
% NeuroStruct out across muliple directories.  This function with gather
% the data out of all of the directories and rebuild it back into the
% NeuroStruct in the same form it would take if we ran the process serially
% instead of in parallel.  

if nargin < 4, MRI_flag = 0; end

% convert char input (bash) to numeric if appropriate
if ischar(chanlist)
   eval(sprintf('chanlist = %s;',chanlist)) 
end
if ischar(MRI_flag)
    eval(sprintf('MRI_flag = %s;',MRI_flag))
end

% get the Tanks name from the path
if strcmp(NeurotankPath(end),filesep)
    NeurotankPath(end) = [];
end
brokenTank = regexp(NeurotankPath,filesep,'split');
tankName = brokenTank{end};

% Some fields are identical across all channels (name, block, and block
% length) while others have data from multiple channels combined into
% single vectors.
for ii = 1:length(chanlist)
    
    curr_chan = num2str(chanlist(ii));
%     load([SavePath filesep curr_chan filesep tankName filesep tankName '_sorted.mat'])
    load(fullfile(SavePath, curr_chan, tankName, [tankName '_sorted.mat']));
    
    if ii == 1
        full_struct = NeuroStruct; % this gets all of the fields that are 
                                   % the same across channels and
                                   % preallocates those that will change
    else
        for jj = 1:length(full_struct) % this fills in all of the data that
                                       % changes across channels.
                        
            full_struct(jj).channels = [full_struct(jj).channels; NeuroStruct(jj).channels];
            full_struct(jj).cells    = [full_struct(jj).cells;    NeuroStruct(jj).cells];
            full_struct(jj).unsorted = [full_struct(jj).unsorted; NeuroStruct(jj).unsorted];
            
        end
    end
end

NeuroStruct = full_struct;

save([SavePath tankName '_sorted.mat'],'NeuroStruct')

if MRI_flag
    add_MRINeuroStruct_spikes2mat(SavePath,[tankName '_sorted.mat'])
end
