% Inclass16

%The folder in this repository contains code implementing a Tracking
%algorithm to match cells (or anything else) between successive frames. 
% It is an implemenation of the algorithm described in this paper: 
%
% Sbalzarini IF, Koumoutsakos P (2005) Feature point tracking and trajectory analysis 
% for video imaging in cell biology. J Struct Biol 151:182?195.
%
%The main function for the code is called MatchFrames.m and it takes three
%arguments: 
% 1. A cell array of data called peaks. Each entry of peaks is data for a
% different time point. Each row in this data should be a different object
% (i.e. a cell) and the columns should be x-coordinate, y-coordinate,
% object area, tracking index, fluorescence intensities (could be multiple
% columns). The tracking index can be initialized to -1 in every row. It will
% be filled in by MatchFrames so that its value gives the row where the
% data on the same cell can be found in the next frame. 
%2. a frame number (frame). The function will fill in the 4th column of the
% array in peaks{frame-1} with the row number of the corresponding cell in
% peaks{frame} as described above.
%3. A single parameter for the matching (L). In the current implementation of the algorithm, 
% the meaning of this parameter is that objects further than L pixels apart will never be matched. 

% Continue working with the nfkb movie you worked with in hw4. 

% Part 1. Use the first 2 frames of the movie. Segment them any way you
% like and fill the peaks cell array as described above so that each of the two cells 
% has 6 column matrix with x,y,area,-1,chan1 intensity, chan 2 intensity
reader = bfGetReader('nfkb_movie1.tif');
z = 1; c = 1; t = 1;
ind = reader.getIndex(z-1, c-1, t-1) +1;
img1_t1 = bfGetPlane(reader, ind);
ind = reader.getIndex(z-1, c-1, t) +1;
img1_t2 = bfGetPlane(reader, ind);
figure;
subplot(1,2,1); imshow(img1_t1, [500 2000]);
subplot(1,2,2); imshow(img1_t2, [500 2000]);
figure; imshowpair(imadjust(img1_t1),imadjust(img1_t2));

reader = bfGetReader('nfkb_movie2.tif');
z = 1; c = 1; t = 1;
ind = reader.getIndex(z-1, c-1, t-1) +1;
img2_t1 = bfGetPlane(reader, ind);
ind = reader.getIndex(z-1, c-1, t) +1;
img2_t2 = bfGetPlane(reader, ind);
figure;
subplot(1,2,1); imshow(img2_t1, [500 2000]);
subplot(1,2,2); imshow(img2_t2, [500 2000]);
figure; imshowpair(imadjust(img2_t1),imadjust(img2_t2));

img1_t1 = img1_t1>1000;
img1_t1 = imclose(img1_t1, strel('disk', 10));
img1_t2 = img1_t2>1000;
img1_t2 = imclose(img1_t2, strel('disk', 10));
img2_t1 = img2_t1>1000;
img2_t1 = imclose(img2_t1, strel('disk', 10));
img2_t2 = img2_t2>1000;
img2_t2 = imclose(img2_t2, strel('disk', 10));
stats= regionprops(img1_t1,'Area');
figure; hist([stats.Area],40);
xlabel('Cell Area', 'FontSize', 24);
ylabel('Frequency', 'FontSize', 24);
% The minarea is around 750
%%
minarea = 750;
mask1_t1 = imfill(img1_t1, 'holes');
mask1_t1 = bwareaopen(mask1_t1, minarea);
mask1_t2 = imfill(img1_t2,'holes');
mask1_t2 = bwareaopen(mask1_t2, minarea);
figure; imshowpair(mask1_t1, mask1_t2);
mask2_t1 = imfill(img2_t1, 'holes');
mask2_t1 = bwareaopen(mask2_t1, minarea);
mask2_t2 = imfill(img2_t2,'holes');
mask2_t2 = bwareaopen(mask2_t2, minarea);
figure; imshowpair(mask2_t1, mask2_t2);

stats1_t1 = regionprops(mask1_t1, img1_t1, 'Centroid', 'Area', 'MeanIntensity');
stats1_t2 = regionprops(mask1_t2, img1_t2, 'Centroid', 'Area', 'MeanIntensity');
stats2_t1 = regionprops(mask2_t1, img2_t1, 'Centroid', 'Area', 'MeanIntensity');
stats2_t2 = regionprops(mask2_t2, img2_t2, 'Centroid', 'Area', 'MeanIntensity');
%%
xy1 = cat(1, stats1_t1.Centroid);
a1 = cat(1, stats1_t1.Area);
mi1 = cat(1,stats1_t1.MeanIntensity);
mi2 = cat(1,stats2_t1.MeanIntensity);
tmp = -1*ones(size(a1));
peaks{1}=[xy1,a1, tmp, mi1, mi2];

xy2 = cat(1, stats1_t2.Centroid);
a2 = cat(1, stats1_t2.Area);
mi21 = cat(1,stats1_t2.MeanIntensity);
mi22 = cat(1,stats2_t2.MeanIntensity);
tmp = -1*ones(size(a2));
peaks{1}=[xy1,a1, tmp, mi21, mi22];
peaks
% need to separate x, y but do not know how to do it. so the dimension does
% not match. And thus the peaks cannot access. 


% Part 2. Run match frames on this peaks array. ensure that it has filled
% the entries in peaks as described above. 


% Part 3. Display the image from the second frame. For each cell that was
% matched, plot its position in frame 2 with a blue square, its position in
% frame 1 with a red star, and connect these two with a green line. 

Figure;
imshow(img1_t1, [500,2000]);
hold on;
plot(peaks{1}(:,1), peaks{1}(:,2), 'r*', 'MarkerSize',24);
plot(peaks{2}(:,1), peaks{2}(:,2), 'cs', 'MarkerSize', 24);
