function [rd_vol]=decimate3D(vol,factors,threshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [rd_vol]=decimate3D(vol,factors)
%   decimate3D takes a 3 dimensional matrix VOL (NxMxW)and reduces it's size by the 
%     factors defined in FACTORS (x,y,z).
%   FACTORS should be an array represented how much smaller each dimension
%     should be.
%   RD_VOL is the reduced volume representing the mean value of 3D area it
%     represents from the original volume.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    thres_flag=0;
else
    thres_flag=1;
end
%%
volsize=size(vol); % get the dimensions of the orginal volume

newsize=volsize./factors; % determine the size of the reduced volume

newshape=[factors(1) newsize(1) factors(2) newsize(2) factors(3) newsize(3)];

Newvol=reshape(vol,newshape);
%%
rd_vol=squeeze(sum(sum(sum(Newvol,5),3),1)); % sum the new dimensions and squeeze them

rd_vol= rd_vol ./ (factors(1)*factors(2)*factors(3));  % get the average for the new volume

%%
if thres_flag  % if a threshold was sent, apply the threshold.
    rd_vol(rd_vol>=threshold)=1; % set everything greater than threshold to 1
    rd_vol(rd_vol<threshold)=0;  % set everything less than threshold to 0
else
    rd_vol=rd_vol;  % if no threshold was sent send back scaled ROI values
end
