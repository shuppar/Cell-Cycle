%% A code to count mRNA in cells. It counts the spots automatically. Cell segmentation is done manually using intensities from autofluorescence or from non-specific binding.

% Copyright (c), Shivnarayan Dhuppar

% In section 3, you might want to change parametes for LMask function depending on the the typical size of the nucleus for the cell line.

close all;
clear all;
clc;


%% Reading images

fmain=sprintf('1');
fileFolder = fullfile('1'); 
dirDAPI = dir(fullfile(fileFolder, '*C0001.tif')); 
dirp53 = dir(fullfile(fileFolder, '*C0002.tif'));
dirP53mRNA = dir(fullfile(fileFolder, '*C0003.tif'));
[DAPI_S, DAPI] = AvgNStack(fileFolder, dirDAPI); % July29_2017
[p53_S, p53] = AvgNStack(fileFolder, dirp53);
[P53_mRNA_S, AvgmRNA] = AvgNStack(fileFolder, dirP53mRNA);  % mRNA is stack of mRNA images.
zs = length(dirDAPI);
[nrows, ncols, nstks] = size(P53_mRNA_S); tot = nrows*ncols;


%% Background Subtraction

if exist('bg.dat', 'file')
    
    bagr = load('bg.dat');
    DB = bagr(1);
    RB = bagr(2);
    PB = bagr(3);
    
else

M = imdilate(im2bw(DAPI, graythresh(DAPI)), strel('disk', 121));
M = ~M;
[~, bgo] = bwlabel(M);
% bgo
bgrnd = regionprops(M, 'Area');
Area = 0;

for i = 1:bgo
    Area = Area + bgrnd(i).Area;
end

M = uint16(M);
DB = max(DAPI_S, [], 3).*(M);
RB = max(P53_mRNA_S, [], 3).*(M);
PB = max(p53_S, [], 3).*(M);

DB = (mean(DB(:))*tot)/(Area); RB = (mean(RB(:))*tot)/(Area); PB = (mean(PB(:))*tot)/(Area);


f = fopen('bg.dat', 'w');  
fprintf(f,'%f\t%f\t%f\n', DB, RB, PB);
% clear dbgo; clear kbgo; clear hbgo; clear M;

end

DAPI = DAPI - DB; p53 = p53 - PB;

for i = 1:nstks   
    P53_mRNA_S(:,:,i) = P53_mRNA_S(:,:,i) - RB;    
end

AvgmRNA = uint16(mean(P53_mRNA_S, 3));

CellSeg = mean(P53_mRNA_S, 3);

%% Nuclear segmentation

[L, num] = NMask(DAPI, 41, 10000, 75000, 0.35, 1.4, 0.35, 1.4);  % NMask(Image, DiskRadius, minArea, maxArea, minCircularity, maxCircularity, minRoundness, maxRoundness)
fprintf('\n\nThe number of cells detected is: %d\n\n', num);
NumLab = regionprops(L, 'Centroid', 'PixelIdxList', 'BoundingBox', 'Area');

figure, imshow(DAPI, []); Figure, imshow(L, []);

%% p53 Annular

N = imdilate(L, strel('disk', 15)); % You might want to change the size of the shell depending on the type of cell line.
N = N - L;

%% Laplacian of Gaussian Filter

forspot = imtophat(P53_mRNA_S, strel('disk', 5));  % Added July29_2017
forspot = double(forspot);
forspot = LOG_filter(forspot);
forspot = forspot/max(forspot(:));

%%

threshGv = zeros(71,1);                      % Size on which the for loop is run.
% B = forspot(:,:, (round(zs/2)));
B = max(forspot, [], 3);
Bprime = im2bw(L, graythresh(L));
Bprime = im2double(~Bprime);
B = B.*Bprime; clear Bprime;

for i = 1:71                               % Number of first maxima to be discarded.
    
    [~, idx] = max(B(:));                       % This finds just the first maximum.
    [x, y] = ind2sub(size(B), idx);
    gtbox = B(shuppar((x-10)):dhuppar(nrows, (x+10)), shuppar((y-10)):dhuppar(ncols, (y+10)));
    temp1 = max(gtbox(:));  % Added July29_2017
    gtbox = gtbox/temp1;
    threshGv(i) = graythresh(gtbox)*temp1; clear temp1; % temp1 added July29_2017;
    B(shuppar((x-10)):dhuppar(nrows, (x+10)), shuppar((y-10)):dhuppar(ncols, (y+10))) = 0;
    
end

threshG = mean(threshGv);
spotcount = forspot > threshG;

Figure, imshow(max(P53_mRNA_S, [],3), []); figure, imshow(max(forspot, [],3));
fprintf('\n\nFrom here starts spot counting. You need to segment cells from the intensity from nonspecific binding of probes.\n\n');

%% removing false spots (area criterion)  % Added July29_2017

[temp1, num1] = bwlabeln(spotcount);
chhotamota = regionprops(temp1, 'Area', 'PixelIdxList');

for i = 1:num1
        if chhotamota(i).Area <= 50 ...   % put in your own criteria for the area... read the starting note.
           && chhotamota(i).Area >= 5       
        else
              spotcount(chhotamota(i).PixelIdxList) = 0;
        end    
end

clear temp1; clear num1;

%% Couting number of spots (modified... Nov 16, 2016).

d = regionprops(L, DAPI, 'Area', 'MeanIntensity');
pn = regionprops(L, p53, 'MeanIntensity');
pc = regionprops(N, p53, 'Area', 'MeanIntensity');

countG = zeros(num,1);  % G maane global.
ManG = zeros(num,1);
DI = zeros(num, 1);
PN = zeros(num,1); PC = zeros(num,1);

%%

for k = 1:num
    
    idx1 = d(k).Area;
    idx2 = d(k).MeanIntensity;
    idx3 = pn(k).MeanIntensity;
    idx4 = pc(k).Area; idx5 = pc(k).MeanIntensity;
    
    DI(k) = idx1*idx2;
    PN(k) = idx3*idx1;
            
        
%% Manual Segmentation

            ManSeg = CellSeg;
            ManSeg(NumLab(k).PixelIdxList) = 0;
    
    
                figure, h_im = imshow(ManSeg, []);
                f = imfreehand;
                BW = createMask(f,h_im);
                BW = im2double(BW);
                cellarea = regionprops(BW, 'Area');
            
%% Counting the mRNAs
                        
                    for j = 1:zs
                   
                     ManTemp(:,:,j) = spotcount(:,:,j).*BW;
            
                    end

                    [~, ManG(k)] = bwlabeln(ManTemp);
            
%% Writing the file

                        f = fopen('Count_compare.dat','a');  
                        fprintf(f,'%f\t%f\t%f\t%f\t%f\n', DI(k), PN(k), idx3/idx5, ManG(k), cellarea(1).Area);
    
                        close all; clear ManTemp; clear ManSeg;
    
end



%% Colouring Cells in different stages of cell cycle in different colours

% if exist('G1_peak.dat', 'file')
%  
% peak1 = load('G1_peak.dat'); 
% peak1 = peak1(:,1);
% CellCycleMask = CellCycleMask(C, L, peak1);
% figure, imshow(DAPI, []);
% hold on
% handle = imshow(Mask1);
% hold off
% set(handle, 'AlphaData', 0.3);
% 
%     hold on                             % To display spot count with the cell cycle mask.
%     for k = 1:num
%     numlab = NumLab(k).Centroid;
%     text(numlab(1), numlab(2), sprintf('%d', count(k)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     
%     end
%     hold off
%     
% else
% 
% end


%% Functions

% %% shuppar function. It makes sure that the cropping region does not go beyond the image size.
% 
% function out = shuppar(x)
% 
% if x > 0
%     out = x;
% else
%     out = 1;
%     
% end
% 
% end
% 
% %% dhuppar function. It makes sure that the cropping region does not go beyond the image size.
% 
% function out = dhuppar(nrows, x)
% 
% if nrows - x > 0
%     out = x;
% else
%     out = nrows;
%     
% end
% end

%%%%%%%%%%% The End %%%%%%%%%%%%%%%%%%%%
