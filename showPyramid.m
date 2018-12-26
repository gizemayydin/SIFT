function showPyramid(numOctaves, numElements, imgPyramid)
plotDim = ceil(sqrt(numElements));
for i=1:1:numOctaves
    figure;
    for j=1:1:numElements
    subplot(plotDim,plotDim,j);
    imshow(imgPyramid{i,j});
    end
end
end