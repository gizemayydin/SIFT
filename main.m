%read the image
img = imread('cameraman.jpg');
%img = rgb2gray(img);

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
          
%show the pyramid
for i=1:1:4
    figure;
    for j=1:1:4
    subplot(2,2,j);
    imshow(dogPyramid{i,j});
    end
end
       
%local max and min
%x,y,value,octave no,row in octav, scale
localMax = [];
localMin = [];

for i=1:1:4 %for octav
    for j=2:1:3 %row in octav
        myImg = dogPyramid{i,j};
        upImg = dogPyramid{i,j-1};
        downImg = dogPyramid{i,j+1};
        rows = size(myImg,1);
        cols = size(myImg,2);
        for k=2:1:rows-1 %row of image
            for l=2:1:cols-1 %col of image
                
                %check neighbors
                value = myImg(k,l);
                neighbors = [myImg(k-1,l-1) myImg(k-1,l) myImg(k-1,l+1) myImg(k,l-1) myImg(k,l+1) myImg(k+1,l-1) myImg(k+1,l) myImg(k+1,l+1) ...
                              upImg(k-1,l-1) upImg(k-1,l) upImg(k-1,l+1) upImg(k,l-1) upImg(k,l) upImg(k,l+1) upImg(k+1,l-1) upImg(k+1,l) upImg(k+1,l+1) ...
                              downImg(k-1,l-1) downImg(k-1,l) downImg(k-1,l+1) downImg(k,l-1) downImg(k,l) downImg(k,l+1) downImg(k+1,l-1) downImg(k+1,l) downImg(k+1,l+1)];
                ismax = 1;
                ismin = 1;
                for m=1:1:26
                    if ismax ==1 || ismin == 1
                        if value <= neighbors(m)
                            ismax = 0;
                        end
                        if value >= neighbors(m)
                            ismin = 0;
                        end
                    else
                        break;
                    end
                end
                
                %it is a max
                if ismax == 1 && ismin == 0
                    localMax = [localMax; k l value i j scales(i,j)];
                end
                %it is a min
                if ismin == 1 && ismax ==0
                    localMin = [localMin; k l value i j scales(i,j)];
                end
            end
        end
    end
end
