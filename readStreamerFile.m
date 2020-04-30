function sev = readStreamerFile(sevFile);
%
% sev = readStreamerFile(sevFile);
%
% Opens a .sev file (raw data saved to TDT streamer) and returns a vector
% of floating point voltage readings (in volts).
%
% last modified 2012-nov-08
% dbtm

fid = fopen(sevFile);
sev = fread(fid,'float32');
fclose(fid);