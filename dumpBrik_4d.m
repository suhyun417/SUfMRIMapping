% dumpBrik_4d.m

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
% addpath(fullfile(dirProjects, 'parksh/_toolbox/Boot_Time_Series'))

% Set directories 
nameSubjNeural = 'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);


pname = [dirDataBOLD, '/tempSURF/'];
dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];



vol = cat(4, dataBOLD.mvoltc{1:3});

% vol = single(DSP.proc.fncvol_3d);


prefix = 'tmp';
filedefault = sprintf('./%s+orig.BRIK',prefix);

[fname,pname] = uiputfile('%.BRIK','Dump Functional BRIK Volume',filedefault);

Info = DATA.Info;
Info.DATASET_RANK = DATA.Info.DATASET_RANK;
Info.DATASET_DIMENSIONS = DATA.Info.DATASET_DIMENSIONS;
Info.DELTA = DATA.Info.DELTA;
Info.BYTEORDER_STRING = DATA.Info.BYTEORDER_STRING;
Info.TYPESTRING = DATA.Info.TYPESTRING;
Info.SCENE_DATA = DATA.Info.SCENE_DATA;
Info.ORIGIN = DATA.Info.ORIGIN;
Info.TAXIS_FLOATS = DATA.Info.TAXIS_FLOATS;
Info.TAXIS_OFFSETS = DATA.Info.TAXIS_OFFSETS;
Info.IDCODE_STRING = DATA.Info.IDCODE_STRING;
Info.IDCODE_DATE = DATA.Info.IDCODE_DATE;
Info.HISTORY_NOTE = DATA.Info.HISTORY_NOTE;


Info.DATASET_DIMENSIONS = [size(vol) 0 0];
Info.DATASET_RANK(2) = size(vol,4); %1;
Info.BRICK_TYPES = 1*ones(1,size(vol,4));  %3 = float
Info.BRICK_FLOAT_FACS = zeros(1,size(vol,4)); %0;
% Info.BRICK_STATS = []; %[min(vol(:)) max(vol(:))]; % can contain information about the min and max values
% Info.BRICK_LABS = DATA.Info.BRICK_LABS(1:2);
% Info.BRICK_KEYWORDS = DATA.Info.BRICK_KEYWORDS;
Info.TAXIS_NUMS(1) = size(vol,4); %1;
% Info.TypeName = 'float';
Info.TypeBytes = 4; %2?
% Info.ByteOrder = 'ieee-le';
Info.Orientation = DATA.Info.Orientation;
Info.FileFormat = DATA.Info.FileFormat;
Info.RootName = DATA.Info.RootName;
Info.Extension_1D = DATA.Info.Extension_1D;


% Info.TAXIS_OFFSETS = []; % If TAXIS_NIMS[1]>0, array of time offsets
%   Info.IDCODE_STRING = []; % will be added if not present
%   Info.IDCODE_DATE   = []; % hardly used
%   Info.BYTEORDER_STRING = 'LSB_FIRST';
%   Info.BRIK_STATS = []; % can contain informatio about the min and max values
%   Info.BRICK_TYPES = 1*ones(1,t); % 1 = short (recommended), 3 = float
%   Info.BRIK_FLOAT_FACS = []; % probably not needed
%   Info.LABEL_1 = filename;
%   Info.DATASET_NAME = sprintf('./%s',filename);
%   Info.TypeName = 'short';
%   Info.TypeBytes = 2;
%   Info.FileFormat = 'BRIK';
  
  
code0 = getAFNIOrientationCode('LR');
code1 = getAFNIOrientationCode('AP');
code2 = getAFNIOrientationCode('DV');

Info.ORIENT_SPECIFIC = [code0 code1 code2];

doti = strfind(fname,'+');
prefix = fname(1:doti-1);
Opt.Prefix = prefix;
Opt.NoCheck = 0;

fprintf('Writing %s as BRIK and HEAD files to %s\n',prefix,pname);
curdir = cd(pname); 

%
% Compensating for an AFNI inconsistency
%
if isfield(Info,'WARP_TYPE'), 
    Info = rmfield(Info,'WARP_TYPE'); 
end

WriteBrik(vol,Info,Opt);  % needs to have a modified version to overwrite
                        % an existing file with "w+"
                        
%% ADDED BY BER 4/5/13 for alignment purposes.   
% copy data to new file so original dimensions are saved.   
copyDataCMD = ['3dcopy  ' prefix ' ' prefix '_refit'];
disp(copyDataCMD);
system(copyDataCMD);

% refit the x y z dimensioned based on strange AFNI conventions 
refitDataCMD = ['3drefit -xdel 1.5 -ydel 1.5 -zdel 1.5 ' prefix '_refit+orig'];
disp(refitDataCMD);
system(refitDataCMD);

monk=DSP.afni_prefix(1:3);
% now that data is refit, align the center with the surface anatomy center
alignDataCMD = ['@Align_Centers -base ' STDPATH.afni_uw '/Anatomy/' monk '_t1_refit_shft+orig. -dset ' prefix '_refit+orig.'];
disp(alignDataCMD);
system(alignDataCMD);

% finally shift the origin slightly because of strange ANFI stuff
shiftDataCMD = ['3drefit -dyorigin -0.5 -dzorigin 0.5 ' prefix '_refit_shft+orig.'];
disp(shiftDataCMD);
system(shiftDataCMD);
%%

 % Copy relevant files to where the surface files are    
 [s,mess,messID]=movefile([pname, '*_refit_shft+orig*'],  dirSPEC);
 if ~s
     fprintf(1,' %s \n', mess);
 end