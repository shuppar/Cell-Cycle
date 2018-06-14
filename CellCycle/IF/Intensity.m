%% A code to stage cells in respective cell cycle stages from the nuclear intensities
% of DNA binding dyes taken as a proxy for DNA content. The same code is used to quantify nuclear 
% content of proteins.

% Copyright (c), Shivnarayan Dhuppar

% In section 5, you might want to change parametes for NMask function depending on the the typical size of the nucleus for the cell line.

%% Reading images

close all;
clear all;
clc;

fmain=sprintf('1');
fileFolder = fullfile('1');
dirDAPI = dir(fullfile(fileFolder, '*C0001.tif')); 
dirp53 = dir(fullfile(fileFolder, '*C0002.tif'));
dirHistxPh = dir(fullfile(fileFolder, '*C0003.tif'));
[DAPI_S, DAPI] = AvgNStack(fileFolder, dirDAPI);
[p53_S, B] = AvgNStack(fileFolder, dirp53);
[hist_S, C] = AvgNStack(fileFolder, dirHistxPh);

[xs, ys, zs] = size(DAPI_S); tot = xs*ys; clear xs; clear ys;
[nrows, ncols] = size(C);

%% Nonuniform illumination correction. Function saved in home/MatlabCode. Copy of Image Analyst's code. (August 2, 2016)

BlankImage = imread('BlankImage.tif');
DAPI = BackgroundCorrect(DAPI, BlankImage);         % Take a blank fluorescein image.
B = BackgroundCorrect(B, BlankImage);
C = BackgroundCorrect(C, BlankImage);

%% Rolling ball subtraction

% DAPI = imtophat(DAPI, strel('disk', 121));
% B = imtophat(B, strel('disk', 121));
% C = imtophat(C, strel('disk', 121));
% A = DAPI;

%% Subtracting Background

if exist('bg.dat', 'file')
    
    bagr = load('bg.dat');
    DB = bagr(1); KB = bagr(2); HB = bagr(3);
    
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
KB = max(p53_S, [], 3).*(M);
HB = max(hist_S, [], 3).*(M);

DB = (mean(DB(:))*tot)/(Area); KB = (mean(KB(:))*tot)/(Area); HB = (mean(HB(:))*tot)/(Area);

f = fopen('bg.dat', 'w');  
fprintf(f,'%f\t%f\t%f\n', DB, KB, HB);

end

for i = 1:zs
    
    hist_S(:,:,i) = hist_S(:,:,i) - HB;
    
end

C = mean(hist_S, 3);

DAPI = DAPI - DB; B = B - KB;
A = DAPI;

%% Nuclear Mask

[L, num] = NMask(A, 29, 7000, 25000, 0.45, 1.45, 0.45, 1.45);  % NMask(Image, DiskRadius, minArea, maxArea, minCircularity, maxCircularity, minRoundness, maxRoundness)

% num = numel(cm1);
fprintf('\n\nThe number of cells detected is: %d\n\n', num);
NumLab = regionprops(L, 'Centroid', 'PixelIdxList', 'BoundingBox', 'Area');
figure, imshow(L, []); figure, imshow(DAPI,[]);
    
%% H2A.X Foci Count

temp = max(hist_S, [], 3);
temp = imgaussfilt(temp, 1);
temp = imtophat(temp, strel('disk', 12));
% LOG filter
LF = -fspecial('log', 10, 1);
temp = imfilter(temp,LF,'replicate');
temp(find(temp<0)) = 0;
% LOG ends
temp = im2double(temp);
%% 
threshGv = zeros(51,1);                      % Size on which the for loop is run.
% B = forspot(:,:, (round(zs/2)));
temp1 = temp;

for i = 1:51                               % Number of first maxima to be discarded.
    
    [~, idx] = max(temp1(:));                       % This finds just the first maximum.
    [x, y] = ind2sub(size(temp1), idx);
    gtbox = temp1(shuppar((x-20)):dhuppar(nrows, (x+20)), shuppar((y-20)):dhuppar(ncols, (y+20)));
    temp2 = max(gtbox(:));  % Added July29_2017
    gtbox = gtbox/temp2;
    threshGv(i) = graythresh(gtbox)*temp2; clear temp2; % temp1 added July29_2017;
    temp1(shuppar((x-20)):dhuppar(nrows, (x+20)), shuppar((y-20)):dhuppar(ncols, (y+20))) = 0;
    
end

threshG = mean(threshGv);
temp = temp > threshG;
% clear temp1;


%%
temp = imfill(temp, 'holes');
temp = bwareaopen(temp, 7);

%% Removing false spots

% [temp1, num1] = bwlabel(temp);
% temp2 = regionprops(temp1, 'Area', 'PixelIdxList', 'Perimeter', 'MajorAxisLength');
% for i = 1:num1
%     if temp2(i).Area <= 150 ...
%         && (temp2(i).Perimeter)^2/(4*pi*temp2(i).Area) <= 1.6 ...
%         && (temp2(i).Perimeter)^2/(4*pi*temp2(i).Area) >= 0.5
%     
%     else
%         temp(temp2(i).PixelIdxList) = 0;
%     end
% end
% clear temp1; clear num1; clear temp2;

%% Visualizing spots

% vs = bwperim(L);
% vs = vs + temp;
% figure, imshow(vs);

%% Overlaying the cells detected over the original image. (Comment this section and the next when running the code for the set of images through shell script.)

% figure(16), imshow(DAPI, []);                             
% hold on
% handle = imshow(visMask);
% hold off
% set(handle, 'AlphaData', 0.3);
% 
% %% Displaying the numbers over the detected nuclei (Changed it to spot numbers, can change it back to nucleus number. Nov2, 2016)
% hold on
% for k = 1:num
%     numlab = NumLab(k).Centroid;
%     text(numlab(1), numlab(2), sprintf('%d', k), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     
% end
% hold off


%% p53 Annular

N = imdilate(L, strel('disk', 15));
N = N - L;    

%% Getting mean intensity of DAPI, H2AX and Ki67 for the detected cells.

d = regionprops(L, DAPI, 'Area', 'MeanIntensity');
pn = regionprops(L, B, 'MeanIntensity');
pc = regionprops(N, B, 'Area', 'MeanIntensity');

h = regionprops(L, C, 'Area','MeanIntensity');



%% Writing Data

    DI = zeros(num, 1);
    Spot_Hp_no = DI;
    HI = DI;
    PN = DI;
    PC = DI;
    numFoci = DI;
    

for j = 1:num
    
    idx1 = d(j).MeanIntensity;
    idx2 = d(j).Area;
    idx3 = pn(j).MeanIntensity;
    idx4 = pc(j).MeanIntensity;
    idx5 = h(j).MeanIntensity;
    
    DI(j) = idx1*idx2;
    PN(j) = idx3*idx2;
    HI(j) = idx5*idx2;  
  
    %% H2AX Foci Count
    
    temp1 = zeros(size(temp));
    temp1(NumLab(j).PixelIdxList) = 1;
    temp1 = temp1.*temp;
    x1 = NumLab(j).BoundingBox;
    x = x1(1); y = x1(2); wi = x1(3); hi = x1(4); clear x1;
    temp1 = temp1(shuppar((y-10)):dhuppar(nrows, (y+hi+10)), shuppar((x-10)):dhuppar(ncols, (x+wi+10)));
    bgmark = -bwdist(~temp1);
    bgmark(~temp1) = Inf;
    bgmark = watershed(bgmark);
    bgmark(~temp1) = 0;
    temp1 = bgmark; clear bgmark;
    [~, numFoci(j)] = bwlabel(temp1); clear temp1;
    
    %%
    
    f = fopen('Intensity.dat','a');  
    fprintf(f,'%f\t%f\t%f\t%f\t%f\n', DI(j), HI(j), numFoci(j), PN(j), idx3/idx4);
    
    clear idx1; clear idx2; clear idx3; clear idx4;
    
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

%% shuppar function. It makes sure that the cropping region does not go beyond the image size.

% function out = shuppar(x)
% 
% if x > 0
%     out = x;
% else
%     out = 1;
%     
% end
% end

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

%%%%%%%%%%%%%%%%% The end! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



