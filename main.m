%read the image
img = imread('house.jpg');
img = rgb2gray(img);

%expand input image by a factor of two using bilinear interpolation
img = imresize(img,2,'bilinear');
img = imgaussfilt(img,1);
imshow(img);

%scales
scales = [0 sqrt(2) 2 2*sqrt(2) 4;
          sqrt(2) 2 2*sqrt(2) 4 4*sqrt(2);
          2 2*sqrt(2) 4 4*sqrt(2) 8;
          2*sqrt(2) 4 4*sqrt(2) 8 8*sqrt(2)];
      
%resize to form images for the octaves
img2 = imresize(img,0.5,'bilinear');
img3 = imresize(img2,0.5,'bilinear');
img4 = imresize(img3,0.5,'bilinear');

%image pyramid
imgPyramid = {img, imgaussfilt(img,scales(1,2)), imgaussfilt(img,scales(1,3)),imgaussfilt(img,scales(1,4)),imgaussfilt(img,scales(1,5));
           imgaussfilt(img2,scales(2,1)), imgaussfilt(img2,scales(2,2)), imgaussfilt(img2,scales(2,3)),imgaussfilt(img2,scales(2,4)),imgaussfilt(img2,scales(2,5));
           imgaussfilt(img3,scales(3,1)),imgaussfilt(img3,scales(3,2)),imgaussfilt(img3,scales(3,3)),imgaussfilt(img3,scales(3,4)),imgaussfilt(img3,scales(3,5));
           imgaussfilt(img4,scales(4,1)),imgaussfilt(img4,scales(4,2)),imgaussfilt(img4,scales(4,3)),imgaussfilt(img4,scales(4,4)),imgaussfilt(img4,scales(4,5))};

%show the pyramid
for i=1:1:4
    figure;
    for j=1:1:5
    subplot(3,3,j);
    imshow(imgPyramid{i,j});
    end
end

%DOG pyramid
dogPyramid = {imgPyramid{1,2}-imgPyramid{1,1}, imgPyramid{1,3}-imgPyramid{1,2}, imgPyramid{1,4}-imgPyramid{1,3},imgPyramid{1,5}-imgPyramid{1,4};
              imgPyramid{2,2}-imgPyramid{2,1}, imgPyramid{2,3}-imgPyramid{2,2}, imgPyramid{2,4}-imgPyramid{2,3},imgPyramid{2,5}-imgPyramid{2,4};
              imgPyramid{3,2}-imgPyramid{3,1}, imgPyramid{3,3}-imgPyramid{3,2}, imgPyramid{3,4}-imgPyramid{3,3},imgPyramid{3,5}-imgPyramid{3,4};
              imgPyramid{4,2}-imgPyramid{4,1}, imgPyramid{4,3}-imgPyramid{4,2}, imgPyramid{4,4}-imgPyramid{4,3},imgPyramid{4,5}-imgPyramid{4,4}};
          
%scales of DOG, bunu dolduralim
scalesDOG = [];
          
%show the pyramid
for i=1:1:4
    figure;
    for j=1:1:4
    subplot(2,2,j);
    imshow(dogPyramid{i,j});
    end
end
       
%local max and min
%x,y,value,octave no, scale
localMax = [];
localMin = [];

for i=1:1:4 %for octav
    for j=2:1:3 %row in octav
        myImg = dogPyramid{i,j};
        rows = size(myImg,1);
        cols = size(myImg,2);
        for k=2:1:row-1 %row of image
            for l=2:1:cols-1 %col of image
            end
        end
    end
end

% %we will use 4 octaves and 5 blur levels.
% A = imgaussfilt(img,sqrt(2));
% B = imgaussfilt(A,sqrt(2));
% octave1 = [img A B imgaussfilt(img,2*sqrt(2)) imgaussfilt(img,4)];
% set1 = {img,A,B,imgaussfilt(img,2*sqrt(2)),imgaussfilt(img,4)};
% 
% B2 = imresize(B,2/3,'bilinear');
% C = imgaussfilt(B2,sqrt(2));
% octave2 = [B2 C imgaussfilt(C,2) imgaussfilt(C,2*sqrt(2)) imgaussfilt(C,4)];
% set2 = {B2,C,imgaussfilt(C,2),imgaussfilt(C,2*sqrt(2)),imgaussfilt(C,4)};
% 
% C2 = imresize(C,2/3,'bilinear');
% D = imgaussfilt(C2, sqrt(2));
% octave3 = [C2 D imgaussfilt(D,2) imgaussfilt(D,2*sqrt(2)) imgaussfilt(D,4)];
% set3 = {C2,D,imgaussfilt(D,2),imgaussfilt(D,2*sqrt(2)),imgaussfilt(D,4)};
% 
% D2 = imresize(D,2/3,'bilinear');
% E = imgaussfilt(D2, sqrt(2));
% octave4 = [D2 E imgaussfilt(E,2) imgaussfilt(E,2*sqrt(2)) imgaussfilt(E,4)];
% set4 = {D2,E,imgaussfilt(E,2) ,imgaussfilt(E,2*sqrt(2)),imgaussfilt(E,4)};
% 
% %Dog filtered versions
% dog1 = [(set1{1} - set1{2}) (set1{2} - set1{3}) (set1{3} - set1{4}) (set1{4} - set1{5})];
% dog2 = [(set2{1} - set2{2}) (set2{2} - set2{3}) (set2{3} - set2{4}) (set2{4} - set2{5})];
% dog3 = [(set3{1} - set3{2}) (set3{2} - set3{3}) (set3{3} - set3{4}) (set3{4} - set3{5})];
% dog4 = [(set4{1} - set4{2}) (set4{2} - set4{3}) (set4{3} - set4{4}) (set4{4} - set4{5})];
% 
% min_dog1 = [];
% max_dog1 = [];
% min_dog1_loc = [];
% max_dog1_loc = [];
% 
% 
% for i=2:1:1439
%     for j=2:1:10239
%         neighbors = dog1(i-1:i+1,j-1:j+1);
%         if dog1(i,j) < dog1(i-1,j-1) && dog1(i,j) < dog1(i-1,j) && dog1(i,j) < dog1(i-1,j+1) && dog1(i,j) < dog1(i,j-1) && dog1(i,j) < dog1(i,j+1) && dog1(i,j) < dog1(i+1,j-1) && dog1(i,j) < dog1(i+1,j) && dog1(i,j) < dog1(i+1,j+1) 
%             min_dog1 = [min_dog1 dog1(i,j)];
%             min_dog1_loc = [min_dog1_loc; i j];
%         elseif dog1(i,j) > dog1(i-1,j-1) && dog1(i,j) > dog1(i-1,j) && dog1(i,j) > dog1(i-1,j+1) && dog1(i,j) > dog1(i,j-1) && dog1(i,j) > dog1(i,j+1) && dog1(i,j) > dog1(i+1,j-1) && dog1(i,j) > dog1(i+1,j) && dog1(i,j) > dog1(i+1,j+1) 
%             max_dog1 = [max_dog1 dog1(i,j)];
%             max_dog1_loc = [max_dog1_loc; i j];
%         end
%     end
% end
% 
% min_search_loc = uint8((2/3)*min_dog1_loc);
% max_search_loc = (2/3)*max_dog1_loc;
% min_dog2 = [];
% max_dog2 = [];
% min_dog2_loc = [];
% max_dog2_loc = [];    
% 
% 
% for k=1:1:4100
%     i = min_search_loc(1,k);
%     i = min_search_loc(2,k);
%     neighbors = dog2(i-1:i+1,j-1:j+1);
%     if dog2(i,j) < dog2(i-1,j-1) && dog2(i,j) < dog2(i-1,j) && dog2(i,j) < dog2(i-1,j+1) && dog2(i,j) < dog2(i,j-1) && dog2(i,j) < dog2(i,j+1) && dog2(i,j) < dog2(i+1,j-1) && dog2(i,j) < dog2(i+1,j) && dog2(i,j) < dog2(i+1,j+1) 
%         min_dog2 = [min_dog2 dog2(i,j)];
%         min_dog2_loc = [min_dog2_loc; i j];
%     end
% end
% 
% for k=1:1:size(max_search_loc,1)
%     [i j] = max_search_loc(k);
%     neighbors = dog2(i-1:i+1,j-1:j+1);
%     if dog2(i,j) > dog2(i-1,j-1) && dog2(i,j) > dog2(i-1,j) && dog2(i,j) > dog2(i-1,j+1) && dog2(i,j) > dog2(i,j-1) && dog2(i,j) > dog2(i,j+1) && dog2(i,j) > dog2(i+1,j-1) && dog2(i,j) > dog2(i+1,j) && dog2(i,j) > dog2(i+1,j+1) 
%         max_dog2 = [max_dog2 dog2(i,j)];
%         max_dog2_loc = [max_dog2_loc; i j];
%     end
% end
% 
% %Finding maxima and minima
% % for i=2:1:1439
% %     for j=2:1:10239
% %         neighbors = dog1(i-1:i+1,j-1:j+1);
% %         if min(neighbors) == dog1(i,j)
% %             min_dog1 = [min_dog1 dog1(i,j)];
% %             min_dog1_loc = [min_dog1_loc; i j];
% %         elseif max(neighbors) == dog1(i,j)
% %             max_dog1 = [max_dog1 dog1(i,j)];
% %             max_dog1_loc = [max_dog1_loc; i j];
% %         end
% %     end
% % end
% 
% 
% 
% 
% 
%         
% 
% 
