function matchSiftFeatures(img1_path, img2_path,descriptors1,points1,descriptors2,points2,mode)
%if mode = 1, the images are montaged, if it is 2 they are shown on top of
%each other

%Read the image, double its size as it mathches the first octave, and
%smooth with sigma = 1
img1 = imread(img1_path);
img1 = imresize(img1,2,'bilinear');
img1 = imgaussfilt(img1,1);

%Do the same for the second image
img2 = imread(img2_path);
img2 = imresize(img2,2,'bilinear');
img2 = imgaussfilt(img2,1);

%The matchFeatures function matches the descriptors of two images, and
%stores tha matched indexes in indexPairs. The second output variable will
%not be used.
[indexPairs,matchmetric] = matchFeatures(descriptors1,descriptors2);

%matchedLocations gets the points with respect to matched indexes of the
%image. The points arrays do not only include the coordinates, but also the
%octaves and scales the points were obtained from.
matchLocations1 = points1(indexPairs(:,1),:);
matchLocations2 = points2(indexPairs(:,2),:);


%since each octave had different sizes of images, every points' coordinates
%were altered as if they were from the first octave, that is the double
%of size of the original image. After that, the columns were interchanged.
loc1 = [];
for i=1:1:size(matchLocations1,1)
    loc1 = [loc1; uint32(matchLocations1(i,1) * (2^(matchLocations1(i,4)-1))) uint32(matchLocations1(i,2) * (2^(matchLocations1(i,4)-1)))];
end
if isempty(loc1)
    disp("No matches");
    return;
end
loc11 = [loc1(:,2) loc1(:,1) ];

loc2 = [];
for j=1:1:size(matchLocations2,1)
    loc2 = [loc2; uint32(matchLocations2(j,1) * (2^(matchLocations2(j,4)-1))) uint32(matchLocations2(j,2) * (2^(matchLocations2(j,4)-1)))];
end
if isempty(loc2)
    disp("No matches");
    return;
end
loc22 = [loc2(:,2) loc2(:,1) ];

%Visualize the results
if mode == 1
figure;showMatchedFeatures(img1,img2,loc11,loc22,'montage');
elseif mode == 2
figure;showMatchedFeatures(img1,img2,loc11,loc22);
end

end
