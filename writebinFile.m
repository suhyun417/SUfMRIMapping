function writebinFile(bin,filename);
%
% writebinFile(bin,filename);
%
% Writes the values in vector BIN to FILENAME. The following code will
% produce a .bin file readable by Plexon OfflineSorter:
%
% sev = readStreamerFile(sevFile);
% bin = streamer2bin(sev);
% writebinFile(bin,binFile);
%
% last modified 2012-nov-08
% dbtm

fid = fopen(filename, 'w');
fwrite(fid,bin,'int16');  % written column by column.  (use bit12 if 12 bits (signed) are required
fclose(fid);