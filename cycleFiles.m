function cycleFiles()


dirData = '/einstein0/USRlab/data/parksh/'; %/rmov/';
W = what(dirData);

cellIDs = {'003a'};
movIDs  = [1 2 3 7 8 9];

for i=1:length(W.mat)-1
  fNameNeural = W.mat{i+1}
  movID = str2num(fNameNeural(5));
  dotindx = strfind(fNameNeural,'.');
  cellID = fNameNeural([dotindx-4:dotindx-1]);
  
  if ~ismember(cellID,cellIDs), continue; end
  if ~ismember(movID,movIDs), continue; end
  
  showRasters(fNameNeural);

end
