function matchSiftFeatures(img1_path, img2_path,descriptors1,points1,descriptors2,points2)
img1 = imread(img1_path);
img1 = imresize(img1,2,'bilinear');
img1 = imgaussfilt(img1,1);

img2 = imread(img2_path);
img2 = imresize(img2,2,'bilinear');
img2 = imgaussfilt(img2,1);

[indexPairs,matchmetric] = matchFeatures(descriptors1,descriptors2);

matchedLoc1 = points1(indexPairs(:,1),:);

loc1 = [];
for i=1:1:size(matchedLoc1,1)
    loc1 = [loc1; uint32(matchedLoc1(i,1) * (2^(matchedLoc1(i,4)-1))) uint32(matchedLoc1(i,2) * (2^(matchedLoc1(i,4)-1)))];
end
loc11 = [loc1(:,2) loc1(:,1) ];
matchedLoc2 = points2(indexPairs(:,2),:);


loc2 = [];
for j=1:1:size(matchedLoc2,1)
    loc2 = [loc2; uint32(matchedLoc2(j,1) * (2^(matchedLoc2(j,4)-1))) uint32(matchedLoc2(j,2) * (2^(matchedLoc2(j,4)-1)))];
end
loc22 = [loc2(:,2) loc2(:,1) ];

figure;showMatchedFeatures(img1,img2,loc11,loc22,'montage');

end
