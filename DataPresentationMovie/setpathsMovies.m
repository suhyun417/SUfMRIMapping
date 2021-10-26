function pp = setpathsMovies_heba(varargin)
%
%  pp = setpathsMovies
%
% Returns paths to data directories for movies project. 
% Runs necessary addpath commands for analysis on helix, einstein,
% mac, or windows PC.
%
% Input argument PROJECT must be a character sting specifying a
% valid project. Use setpaths('h') to see a list of valid options.
% 
% Output structure fields (and example)
% pp.dataHome   /data/mcmahond/movies
% pp.raw        /rawdata/mcmahond/
% pp.rare       /data/mcmahond/movies/rare
% pp.mas        /data/mcmahond/movies/mas
% pp.ofs        ~/ana/OfflineSDK
% pp.han        ~/ana/han06
% pp.ana        ~/ana/movies
% pp.results    ~/results/movies
% pp.stim       ~/sim/movies
%
% Additional field may be included on einstein and spike sorting PCs:
% pp.preproc    /einstein0/USRlab/projects/mcmahond/preproc/
%
% pp = setpaths(project,OPTION)
%
% OPTION arguement must be a character string. Valid options may be
% used in combination and include:
% 'v' verbose
% 'a' executes matlab addpath function for a subset of paths
% 'r' executes matlab rmpath function for the same subset of paths.
% 
% last modified 2013-july-19
% hde

proj = 'movies';

if isempty(varargin)
    verbose = 0;
    addFlag = 0;
    rmFlag = 0;
else
    option = varargin{1};
    if ~isempty(strfind(option,'v'))
        verbose = 1;
    else
        verbose = 0;
    end
    if ~isempty(strfind(option,'a'))
        addFlag = 1;
    else
        addFlag = 0;
    end
    if ~isempty(strfind(option,'r'))
        rmFlag = 1;
    else
        rmFlag = 0;
    end
end

if ispc
    pp.project = proj;
    pp.dataHome = ['\\nifstorage1.nimh.nih.gov\USRlab\data\mcmahond\' proj '\'];
    pp.raw      = ['\\nifstorage1.nimh.nih.gov\USRlab\rawdata\mcmahond\'];
    pp.rare     = [pp.dataHome 'rare\'];
    pp.mas      = [pp.dataHome 'mas\'];
    pp.ofs      = ['\\nifstorage1.nimh.nih.gov\USRlab\projects\mcmahond\ana\OfflineSDK\'];
    pp.han      = ['\\nifstorage1.nimh.nih.gov\USRlab\projects\mcmahond\ana\han06\'];
    pp.preproc  = ['\\nifstorage1.nimh.nih.gov\USRlab\projects\mcmahond\ana\preproc\'];
    %pp.anaHome  = ['\\nifstorage1.nimh.nih.gov\USRlab\projects\mcmahond\ana\' proj '\'];
    pp.analysis = ['\\nifstorage1.nimh.nih.gov\USRlab\projects\heba\analysis\' proj '\'];
    pp.results  = ['\\nifstorage1.nimh.nih.gov\USRlab\projects\heba\results\' proj '\'];
    pp.stim  = ['\\nifstorage1.nimh.nih.gov\USRlab\projects\heba\stimuli\' proj '\'];
    pp.dbtmShared = ['\\nifstorage1.nimh.nih.gov\usrlab\projects\mcmahond\heba\'];
% if ispc
%     pp.project = proj;
%     pp.dataHome = ['z:\USRlab\data\mcmahond\' proj '\'];
%     pp.raw      = ['z:\USRlab\rawdata\mcmahond\'];
%     pp.rare     = [pp.dataHome 'rare\'];
%     pp.mas      = [pp.dataHome 'mas\'];
%     pp.ofs      = ['z:\USRlab\projects\mcmahond\ana\OfflineSDK\'];
%     pp.han      = ['z:\USRlab\projects\mcmahond\ana\han06\'];
%     pp.preproc  = ['z:\USRlab\projects\mcmahond\ana\preproc\'];
% %     pp.anaHome  = ['z:\USRlab\projects\mcmahond\ana\' proj '\'];
%     pp.analysis = ['z:\USRlab\projects\heba\analysis\' proj '\'];
%     pp.results  = ['z:\USRlab\projects\heba\results\' proj '\'];
%     pp.stim  = ['z:\USRlab\projects\heba\stimuli\' proj '\'];
elseif isunix
    [junk hostname] = system('hostname');
    %    if isequal(hostname,'helix.nih.gov');
%         if isdir('/data/mcmahond/') 
%             pp.project = proj;
%             pp.dataHome = ['/data/mcmahond/' proj '/'];
%             pp.raw      = [pp.dataHome 'raw/'];
%             pp.rare     = [pp.dataHome 'rare/'];
%             pp.mas      = [pp.dataHome 'mas/'];
%             pp.ofs      = ['~/ana/OfflineSDK/'];
%             pp.han      = ['~/ana/han06/'];
%             pp.anaHome      = ['~/ana/stimscreen/'];            
%             pp.results  = ['/projects/mcmahond/results/' proj '/'];
%             pp.stim  = ['~/stim/' proj '/'];            
%         elseif isdir('/einstein0/')
          if isdir('/einstein0/')
            pp.project = proj;
            pp.dataHome = ['/einstein0/USRlab/data/mcmahond/' proj '/'];
            pp.raw      = [];%[pp.dataHome 'raw/'];
            pp.rare     = [pp.dataHome 'rare/'];
            pp.mas      = [pp.dataHome 'mas/'];
            pp.ofs      = ['/einstein0/USRlab/projects/mcmahond/ana/OfflineSDK/'];
            pp.han      = ['/einstein0/USRlab/projects/mcmahond/ana/han06/'];
%             pp.anaHome      = ['/einstein0/USRlab/projects/mcmahond/ana/' proj '/']; 
            pp.analysis  = ['/einstein0/USRlab/projects/heba/analysis/' proj '/']; 
            pp.results  = ['/einstein0/USRlab/projects/heba/results/' proj '/'];
            pp.stim  = ['/einstein0/USRlab/projects/heba/stimuli/' proj '/'];
            pp.dbtmShared = ['/einstein0/USRlab/projects/mcmahond/heba/'];
            
        end
        %end
elseif ismac
    ;
end

% if addFlag==1
    if (~isdeployed)
    addpath('/einstein0/USRlab/projects/mcmahond/ana/han06/');
    addpath(['/einstein0/USRlab/data/mcmahond/movies/mas/']);
    %addpath(pp.anaHome);
    addpath(['/einstein0/USRlab/projects/heba/analysis/movies/']);
    addpath(['/einstein0/USRlab/projects/heba/stimuli/movies/']);
    addpath('/einstein0/USRlab/projects/mcmahond/ana/OfflineSDK/');
    addpath('/einstein0/USRlab/projects/mcmahond/heba/');
    addpath('/einstein0/USRlab/projects/heba/analysis/cbrewer/');
    end
% elseif rmFlag==1
%     ;
%     %    rmpath(pp.mas);
%     %rmpath(pp.anaHome);
%     %rmpath(pp.stim);
% end

if verbose
    disp(pp)
end
