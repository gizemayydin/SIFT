function showPyramid(numOctaves, numElements, imgPyramid)
%Takes an image pyramid and its dimensions as an input and shows each octave
%in a seperate figure
plotDim = ceil(sqrt(numElements));
for i=1:1:numOctaves
    figure;
    for j=1:1:numElements
        subplot(plotDim,plotDim,j);
        imshow(imgPyramid{i,j});
    end
end
end