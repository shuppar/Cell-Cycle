
function [L, num] = NMask(Image, DiskRadius, minArea, maxArea, minCircularity, maxCircularity, minRoundness, maxRoundness)

%%
Q = anisodiff(Image, 7, 20, 0.25, 2);              % anisodiff(im, niter, kappa, lambda, option)
Q = anisodiff(Image, 7, 20, 0.25, 1);              % anisodiff(im, niter, kappa, lambda, option)
%%
A = uint16(Q);
%%
% A = adapthisteq(A);                                 % Some local a djustments, to be able to detect the dimmer cells. Commented out on March31, 2017.

%%
conli = stretchlim(A);
A = imadjust(A, conli, [0 1], 1.3);          % April 1, 2017. Gamma Correction.
%%
% A = imclearborder(A);                               % Eliminating the objects on the borders. Commented: April 8, 2017.
%%
A = wiener2(A, [3*DiskRadius, 3*DiskRadius]);                          % Removing pixels smaller than the given size.

% figure(1), imshow(A), title('original DAPI image');
%figure(2), imshow(B), title('Original H2AX image');
%figure(3), imshow(C), title('Original Ki67 image');


%% Removing the problem of oversegmentation.

se = strel('disk', DiskRadius);
Ae = imerode(A, se);
Aobr = imreconstruct(Ae, A);
Aobrd = imdilate(Aobr, se);
Aobrcbr = imreconstruct(imcomplement(Aobrd), imcomplement(Aobr));
Aobrcbr = imcomplement(Aobrcbr);



%% Now this will be nothing like Pedro's blog.

hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(imclearborder(Image)), hy, 'replicate');
Ix = imfilter(double(imclearborder(Image)), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);


fgm = imregionalmax(Aobrcbr);
fgm = imclose(fgm, strel('disk', 4));
fgm = imerode(fgm, strel('disk', 4));
fgm = bwareaopen(fgm, 100);
bw = im2bw(Aobrcbr, graythresh(A));
Dist = bwdist(bw);
DL = watershed(Dist);
bgm = DL == 0;
gradmag2 = imimposemin(gradmag, bgm | fgm);
L = watershed(gradmag2);


%% Counting cells and removing under- and over segmented nuclei.

[L, num1] = bwlabel(L);
chhotamota = regionprops(L, 'Area', 'PixelIdxList', 'Perimeter', 'MajorAxisLength');


%% Detecting nuclei with optimum area.

for i = 1:num1
%         fprintf('%d\n', chhotamota(i).Area);
        if chhotamota(i).Area <= maxArea ...                                                          % put in your own criteria for the area... read the starting note.
           && chhotamota(i).Area >= minArea ...    
           && (4*pi*chhotamota(i).Area)/(chhotamota(i).Perimeter)^2 <= maxCircularity ...             % Circularity condition
           && (4*pi*chhotamota(i).Area)/(chhotamota(i).Perimeter)^2 >= minCircularity ...
           && (4*chhotamota(i).Area)/(pi*(chhotamota(i).MajorAxisLength)^2) <= maxRoundness ...       % Roundness.          
           && (4*chhotamota(i).Area)/(pi*(chhotamota(i).MajorAxisLength)^2) >= minRoundness;
%            cm1(i) = i;  %commented on July 6, 2017. Just realised I am not using it anywhere.
           
        else

              L(chhotamota(i).PixelIdxList) = 0;
        end    
end

L = imopen(L, strel('disk', 20));       %April 7, 2017; to smooth out the edges.

L = imclearborder(L, 8);
[L, num] = bwlabel(L);


%% dhuppar function

function out = dhuppar(nrows, x)

if nrows - x > 0
    out = x;
else
    out = nrows;
    
end

%% shuppar function

function out = shuppar(x)

if x > 0
    out = x;
else
    out = 1;
    
end


%% Aniso Diff function


% ANISODIFF - Anisotropic diffusion.
%
% Usage:
%  diff = anisodiff(im, niter, kappa, lambda, option)
%
% Arguments:
%         im     - input image
%         niter  - number of iterations.
%         kappa  - conduction coefficient 20-100 ?
%         lambda - max value of .25 for stability
%         option - 1 Perona Malik diffusion equation No 1
%                  2 Perona Malik diffusion equation No 2
%
% Returns:
%         diff   - diffused image.
%
% kappa controls conduction as a function of gradient.  If kappa is low
% small intensity gradients are able to block conduction and hence diffusion
% across step edges.  A large value reduces the influence of intensity
% gradients on conduction.
%
% lambda controls speed of diffusion (you usually want it at a maximum of
% 0.25)
%
% Diffusion equation 1 favours high contrast edges over low contrast ones.
% Diffusion equation 2 favours wide regions over smaller ones.

% Reference: 
% P. Perona and J. Malik. 
% Scale-space and edge detection using ansotropic diffusion.
% IEEE Transactions on Pattern Analysis and Machine Intelligence, 
% 12(7):629-639, July 1990.
%
% Peter Kovesi  
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk @ csse uwa edu au
% http://www.csse.uwa.edu.au
%
% June 2000  original version.       
% March 2002 corrected diffusion eqn No 2.

function diff = anisodiff(im, niter, kappa, lambda, option)

if ndims(im)==3
  error('Anisodiff only operates on 2D grey-scale images');
end

im = double(im);
[rows,cols] = size(im);
diff = im;
  
for i = 1:niter
%  fprintf('\rIteration %d',i);

  % Construct diffl which is the same as diff but
  % has an extra padding of zeros around it.
  diffl = zeros(rows+2, cols+2);
  diffl(2:rows+1, 2:cols+1) = diff;

  % North, South, East and West differences
  deltaN = diffl(1:rows,2:cols+1)   - diff;
  deltaS = diffl(3:rows+2,2:cols+1) - diff;
  deltaE = diffl(2:rows+1,3:cols+2) - diff;
  deltaW = diffl(2:rows+1,1:cols)   - diff;

  % Conduction

  if option == 1
    cN = exp(-(deltaN/kappa).^2);
    cS = exp(-(deltaS/kappa).^2);
    cE = exp(-(deltaE/kappa).^2);
    cW = exp(-(deltaW/kappa).^2);
  elseif option == 2
    cN = 1./(1 + (deltaN/kappa).^2);
    cS = 1./(1 + (deltaS/kappa).^2);
    cE = 1./(1 + (deltaE/kappa).^2);
    cW = 1./(1 + (deltaW/kappa).^2);
  end

  diff = diff + lambda*(cN.*deltaN + cS.*deltaS + cE.*deltaE + cW.*deltaW);

%  Uncomment the following to see a progression of images
%  subplot(ceil(sqrt(niter)),ceil(sqrt(niter)), i)
%  imagesc(diff), colormap(gray), axis image

end
%fprintf('\n');