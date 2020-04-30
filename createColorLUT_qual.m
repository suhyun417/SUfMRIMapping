% createColorLUT_qual.m
% 2017/01/25 SHP
% Create a color look up table, using qualitative LUTs from colorbrewer2.org
% This is for different ROIs (or parcellations), so here I merged two CLUTs, each for 12 classes

clear all;

% % one example
% pastel = [141,211,199;...
% 255,255,179;...
% 190,186,218;...
% 251,128,114;...
% 128,177,211;...
% 253,180,98;...
% 179,222,105;...
% 252,205,229;...
% % 217,217,217;... % this color is too whitish
% 188,128,189;... 
% 204,235,197;... % mint reappears
% 255,255,51;... %255,237,111;...% yellow reappears, so changed it to darker yellow
% ];
% 
% dark = [166,206,227;...
% 31,120,180;...
% 178,223,138;...
% 51,160,44;...
% 251,154,153;...
% 227,26,28;...
% 253,191,111;...
% 255,127,0;...
% 202,178,214;...
% 106,61,154;...
% % 255,255,153;...
% 177,89,40];
% 
% clut = cat(1, pastel, dark)./255; % should be between 0 and 1

% one example
mat1 = [228	26	28;
55	126	184;
77	175	74;
152	78	163;
255	127	0;
255 217 47; % dark yellow %255	255	255; %white %255	255	51; % yellow was too similar to another yellow in mat2
166	86	40;
247	129	191;
];

mat2 = [141	211	199
255	255	179
190	186	218
251	128	114
128	177	211
253	180	98
179	222	105
252	205	229];

clut = cat(1, mat1, mat2)./255; % should be between 0 and 1

% to use it in AFNI
clut = cat(1, [1 1 1], clut);
ind = (0:1:length(clut)-1)'; 
clut = cat(2, clut, ind);

% write to a text file
% for iK=2:5
%     
dlmwrite(sprintf('colorLUT_qual%d_plusWhite.txt', length(clut)-1), clut)
movefile(sprintf('colorLUT_qual%d_plusWhite.txt', length(clut)-1), '/Volumes/PROCDATA/parksh/Art/Anatomy/_suma')

