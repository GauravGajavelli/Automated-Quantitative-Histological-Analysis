function Function_Final_Program_06112020(Filepath, Segments)

disp('Filepath 2')
disp(Filepath)
disp('Segments #2')
disp(Segments)

cd(Filepath); %points to directory where script is located
myFolder = cd;
filePattern = fullfile(myFolder,'*.jpg'); %shows file path to high resolution images
disp(filePattern);
theFiles = dir(filePattern); %creates an n x 1 array where n is number of rows = # of *.jpg files in the originalhigh res folder
[numRows,numCols] = size(theFiles);

disp(numRows);
%setting the directory to the given source folder 'OriginalResImage'

source = cd; %type in the name of the source folder containing the images

source = [source, '\'];
mkdir(source,'CopyOfImages'); %makes a subfolder to copy the read images
mkdir(source,'Montages'); %makes a subfolder to copy the read images
mkdir(source,'Segmented'); %makes a subfolder to copy the read images

% output = inputdlg('Number of Segments (input of 2 is accepted) (multiples of 4):');
% number = str2num(output{1,1});
segments = Segments;

cd(source);
folders = dir;
disp(folders);


         for k = 1 : length(theFiles) %loop through all images in the given source folder
            
             baseFileName = theFiles(k).name;
             fullFileName = fullfile(baseFileName); %get the name of the current image
             fprintf(1, 'Now reading %s\n', fullFileName);
             I = imread(fullFileName);
             stringName = strcat(fullFileName);
             baseFileName = sprintf(stringName);
             fullFileName2 = fullfile('CopyOfImages', baseFileName);

            %find the x value with the lowest y value in the range of the 2nd and 3rd quartiles
            level = graythresh(I);
            BW = imbinarize(I,level); %Makes image black and white
            imshow(BW(:,:,3:3));
            BW3 = BW(:,:,3:3);
            BW3comp = imcomplement(BW3); %makes white black and black white in the image
            [rows, columns, numberOfColorChannels] = size(BW); %Takes the size of the image and puts in an array of variables
            colsums = zeros(columns,1);
            initial = 1;
            colsum = 0;
            for i = 1:columns-1 %i is for columns
                for j = 1:rows-1 %j is for rows
                    colsum = colsum + BW3comp(j,i); %adds up all the pixels in each column
                end
                colsums(i,1) = colsum;
                %         colsums(i,2) = colsum;
                %         colsums(i,1) = initial;
                %         initial = initial + 1;
                colsum = 0;
            end
           
            
            
            nzero = find(colsums); %puts all non-zero values in an array
            %disp(k)
            first_ele=nzero(1,:);   %first value of nzero
            last_ele=nzero(end,:);   %last value of nzero
            fprintf('The first element: %d\n', first_ele);
            fprintf('The last element: %d\n', last_ele);
            range = last_ele-first_ele; %finds the range of the non-zero elements
            interquart = range/segments;
            
            fprintf('The Range: %d\n', range);
            fprintf('The interquartile Range in eighths: %d\n', interquart);
            
            q1 = round(interquart +first_ele);
            fprintf('Quartile 1: %d\n', q1);
            q2 = round(interquart*2 +first_ele);
            fprintf('Quartile 2: %d\n', q2);
            q3 = round(interquart *3 +first_ele);
            fprintf('Quartile 3: %d\n', q3);
            q4 = round(interquart*4 +first_ele); %gets the second quartile value
            fprintf('Quartile 4: %d\n', q4);
            
           
            q5 = round(interquart*5+first_ele +100); %gets the third quartile value
            fprintf('Quartile 5: %d\n', q5);
            
            q6= round(interquart*6 +first_ele);
            fprintf('Quartile 6: %d\n', q6);
            q7 = round(interquart*7 +first_ele);
            fprintf('Quartile 7: %d\n', q7);
            q8 = round(interquart*8 +first_ele);
            fprintf('Quartile 8: %d\n', q8);
            
            q9 = round(interquart*9 +first_ele);
            fprintf('Quartile 9: %d\n', q9);
            q10 = round(interquart*10 +first_ele);
            fprintf('Quartile 10: %d\n', q10);
            
            if(segments ==2)
            x=0;
            min = colsums(q1); %creates a value for comparison that can be useed to find the column with the lowest y value
            for i = q1 : q2 %i goes through each value in between the quartile values
                if(colsums(i-1)<min) %find the lowest y value
                 min = colsums(i-1);
                    x=i-1;
                end
                i=i+1;
            end
            end
            
            if(segments ==8)
            x=0;
            min = colsums(q3); %creates a value for comparison that can be useed to find the column with the lowest y value
            for i = q3 : q5 %i goes through each value in between the quartile values
                if(colsums(i-1)<min) %find the lowest y value
                 min = colsums(i-1);
                    x=i-1;
                end
                i=i+1;
            end
            end
            if(segments == 4)
            x=0;
            min = colsums(q2); %creates a value for comparison that can be useed to find the column with the lowest y value
            for i = q2 : q3 %i goes through each value in between the quartile values
                if(colsums(i-1)<min) %find the lowest y value
                 min = colsums(i-1);
                    x=i-1;
                end
                i=i+1;
            end
            end
            if(segments == 16)
            x=0;
            min = colsums(q6); %creates a value for comparison that can be useed to find the column with the lowest y value
            for i = q7 : q10 %i goes through each value in between the quartile values
                if(colsums(i-1)<min) %find the lowest y value
                 min = colsums(i-1);
                    x=i-1;
                end
                i=i+1;
            end
            end
            
            
            
            minX = x;
            fprintf('The minimum x value: %d\n', minX);
            fprintf('The minimum y value: %d\n', colsums(minX));
            
            Img = imread(fullFileName); %read the image in again
            

            imgGray=rgb2gray(Img);
            imgNew = rgb2gray(Img); %convert the color image to grayscale
            
            
            imgGray(:, minX) = 0; % White = 255, can pick any intensity. Burns in Line
            figure(1)
            
            if (segments==2)
            t = sprintf('Segment 1 =%f pixel',q1);
            imgGray = insertText(imgGray, [q1,1600], t, 'FontSize', 32);
             imgGray(:, q1) = 0;
            t = sprintf('Segment 2 =%f pixel',q2);
            imgGray = insertText(imgGray, [q2,1700], t, 'FontSize', 32);
             imgGray(:, q2) = 0;
            end
            if (segments >2)
                t = sprintf('Segment 1 =%f pixel',q1);
            imgGray = insertText(imgGray, [q1,1600], t, 'FontSize', 32);
             imgGray(:, q1) = 0;
            t = sprintf('Segment 2 =%f pixel',q2);
            imgGray = insertText(imgGray, [q2,1700], t, 'FontSize', 32);
             imgGray(:, q2) = 0;
            t = sprintf('Segment 3 =%f pixel',q3);
            imgGray = insertText(imgGray, [q3,1800], t, 'FontSize', 32);
             imgGray(:, q3) = 0;
            t = sprintf('Segment 6 =%f pixel',q6);
            imgGray = insertText(imgGray, [q6,2000], t, 'FontSize', 32);
             imgGray(:, q6) = 0;
            t = sprintf('Segment 7 =%f pixel',q7);
            imgGray = insertText(imgGray, [q7,2100], t, 'FontSize', 32);
             imgGray(:, q7) = 0;
            t = sprintf('Segment 8 =%f pixel',q8);
            imgGray = insertText(imgGray, [q8,2300], t, 'FontSize', 32);
             imgGray(:, q8) = 0;
            
            t = sprintf('Segment 4 =%f pixel',q4);
            imgGray = insertText(imgGray, [q4,2200], t, 'FontSize', 32);
             imgGray(:, q4) = 0;
            
            
            t = sprintf('Segment 5 =%f pixel',q5);
            imgGray = insertText(imgGray, [q5,1500], t, 'FontSize', 32);
             imgGray(:, q5) = 0;
            
            end
            
            t = sprintf('The line is at x =%f pixel',minX);
            imgGray = insertText(imgGray, [minX,1900], t, 'FontSize', 32);
            imshow(imgGray);
            h = gca;
            h.Visible = 'On';
            [rows, columns, numberOfColorChannels] = size(Img);
            L = line([minX, minX], [1, rows],'LineStyle','-','Color','k'); %draw in a line to show where we are going to cut the image
            set(L, 'Linewidth', 1);
            
            stringName = strcat(fullFileName, '_grayscale_wLine.jpeg');
            baseFileName = sprintf(stringName);
            fullFileName2 = fullfile('CopyOfImages', baseFileName);

            
            
            
            imwrite(imgGray,fullFileName2)%saving of the grayscale image
            
            caption = sprintf('This is where the lowest concentration of pixels and center are. \n Location = %d pixels',minX);
            title(caption, 'FontSize', 11);
            
            [rows, columns] = size(imgGray);
            middleCol = minX;
            leftHalf = imgNew(:, 1:middleCol); %splits the image into its left half and makes a new image
            rightHalf = imgNew(:, middleCol+1:end); %splits the image into its right half and makes a new image
            
            stringName = strcat(fullFileName, '_leftHalf.jpeg');
            baseFileName = sprintf(stringName);
            fullFileName2 = fullfile('Segmented', baseFileName);
            
            imwrite(leftHalf,'leftHalf_Brain.jpeg', 'quality', 90)
            imwrite(leftHalf, fullFileName2);
            
            stringName = strcat(fullFileName, '_rightHalf.jpeg');
            baseFileName = sprintf(stringName);
            fullFileName2 = fullfile('Segmented', baseFileName);
           
            imwrite(rightHalf,'rightHalf_Brain.jpeg', 'quality', 90)
            imwrite(rightHalf, fullFileName2);
            
            
            
            figure(2)
            chart=bar(colsums); %create the bar chart with each column value of the original image 
           
            %x = size(baseFileName);
            %s = source+"bars\"+baseFileName(1:x(2)-4);
            %writematrix(colsums,s+"Matrix.xlsx")
            
            stringName = strcat(fullFileName, '_BarGraph.jpeg');
            baseFileName = sprintf(stringName);
            fullFileName2 = fullfile('CopyOfImages', baseFileName);
            
            saveas(chart,fullFileName2);
            
            saveas(chart, 'Chart_MontageUse.jpeg');
            ImgChart = imread('Chart_MontageUse.jpeg');
            
            
            
            
            
            
%Begin Analyzation of Left Half

ImgLeft = imread('leftHalf_Brain.jpeg');
figure(3);
imshow(ImgLeft);

%Converts image to Black and white

img2=im2bw(ImgLeft,graythresh(ImgLeft));
imshow(img2)
%Changes black to white and vice versa in image
img2=~img2;
imshow(img2)
%Example 1: Text insertion
%imgOut = insertInImage(im, @()text(40,50,'This is embedded text!'),...
    %{'fontweight','bold','color','m','fontsize',11,...
    %'linewidth',3,'margin',5,'edgecolor',[0 1 1],'backgroundcolor','y'});
%imshow(imgOut);
%imwrite(imgOut, 'textinsertionimage.jpg','quality',90) %the image with the text that is a test
 

%imgOut = insertInImage(imbw, @()text(40,50,'This is embedded text!'),...
    %{'fontweight','bold','color','m','fontsize',11,...
    %'linewidth',3,'margin',5,'edgecolor',[0 1 1],'backgroundcolor','y'});
%imshow(imgOut);
%imwrite(imgOut, 'imgbw2.jpg','quality',90) %the black and white photo with text


%Finds the boundaries for each object and stores in variable B
%Number of Objects found is printed on Image by Text Function  
B = bwboundaries(img2);
imshow(img2)
text(4,4,strcat('\color{green}Objects Found:',num2str(length(B))))
hold on 

%Marks the boundaries of objects with green lines
for cnt = 1:length(B)
boundary = B{cnt};
plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1)
end

imshow(img2)

% Identify individual blobs by seeing which pixels are connected to each other.
% Each group of connected pixels will be given a label, a number, to identify it and distinguish it from the other blobs.
% Do connected components labeling with either bwlabel() or bwconncomp().
labelednImage = bwlabel(img2, 8);     % Label each blob so we can make measurements of it
% labeledImage is an integer-valued image where all pixels in the blobs have values of 1, or 2, or 3, or ... etc.
subplot(3, 3, 1);
imshow(labelednImage, []);  % Show the gray scale image.

%imwrite(labelednImage,'labeledimage1.jpg','quality', 90)
%imlabeled1 = imread('labeledimage1.jpg');
title('Labeled Image, from bwlabel()');

% Let's assign each blob a different color to visually show the user the distinct blobs.
coloredLabels = label2rgb (labelednImage, 'hsv', 'k', 'shuffle'); % pseudo random color labels
% coloredLabels is an RGB image.  We could have applied a colormap instead (but only with R2014b and later)
subplot(3, 3, 5);
imshow(coloredLabels);

axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
caption = sprintf('Pseudo colored labels, from label2rgb().\nBlobs are numbered from top to bottom, then from left to right.');
title(caption);

% Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
blobMeasurements = regionprops(labelednImage, img2, 'all');
numberOfBlobs = size(blobMeasurements, 1);

% bwboundaries() returns a cell array, where each cell contains the row/column coordinates for an object in the image.
% Plot the borders of all the coins on the original grayscale image using the coordinates returned by bwboundaries.
subplot(3, 3, 6);
imshow(img2);
title('Outlines, from bwboundaries()'); 
axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
hold on;
boundaries = bwboundaries(img2,4);
numberOfBoundaries = size(boundaries, 1);
for cnt = 1 : numberOfBoundaries
	thisBoundary = boundaries{cnt};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
end
hold off;

textFontSize = 14;	% Used to control size of "blob number" labels put atop the image.
labelShiftX = -7;	% Used to align the labels in the centers of the coins.
blobECD = zeros(1, numberOfBlobs);

filename = 'Export_GajavelliData.xlsx';% create file to export data
TableHeader = {'Healthy Brain Half','Mean Intensity', 'Area', 'Perimeter', 'Diameter'};
writecell(TableHeader,filename,'Sheet','Data','Range','A1');


% Loop over all blobs printing their measurements to the command window.
for count = 1 : numberOfBlobs           % Loop through all blobs.
	% Find the mean of each blob.  (R2008a has a better way where you can pass the original image
	% directly into regionprops.  The way below works for all versions including earlier versions.)
	thisBlobsPixels = blobMeasurements(count).PixelIdxList;  % Get list of pixels in current blob.
	meanGL = mean(img2(thisBlobsPixels)); % Find mean intensity (in original image!)
	meanGL2008a = blobMeasurements(count).MeanIntensity; % Mean again, but only for version >= R2008a
	
	blobArea = blobMeasurements(count).Area;		% Get area.
	blobPerimeter = blobMeasurements(count).Perimeter;		% Get perimeter.
	blobCentroid = blobMeasurements(count).Centroid;		% Get centroid one at a time
	blobECD(count) = sqrt(4 * blobArea / pi);					% Compute ECD - Equivalent Circular Diameter.
    if blobArea > 5000
    fprintf('The Area of the Left Half: %d\n', blobArea);
    healthymeasurements(1,:) = [1 meanGL blobArea blobPerimeter blobECD(count)]; %Matrix to store collected data
    end
	% Put the "blob number" labels on the "boundaries" grayscale image.
	text(blobCentroid(1) + labelShiftX, blobCentroid(2), num2str(count));
end


%{ 
    Put the labels on the rgb labeled image also.
subplot(3, 3, 5);
for k = 1 : numberOfBlobs           % Loop through all blobs.
	text(centroidsX(k) + labelShiftX, centroidsY(k), num2str(k), 'FontSize', textFontSize, 'FontWeight', 'Bold');
end
%}


allBlobAreas = [blobMeasurements.Area];
% Get a list of the blobs that meet our criteria and we need to keep.
allowableAreaIndexes = allBlobAreas > 2000; % Take the small objects.
keeperIndexes = find(allowableAreaIndexes);
% Extract only those blobs that meet our criteria, and
% eliminate those blobs that don't meet our criteria.
% Note we use ismember() to do this.  Result will be an image - the same as labeledImage but with only the blobs listed in keeperIndexes in it.
keeperBlobsImage = ismember(labelednImage, keeperIndexes);
% Re-label with only the keeper blobs kept.
labelednImage = bwlabel(keeperBlobsImage, 8);     % Label each blob so we can make measurements of it
subplot(3, 3, 7);
imshow(labelednImage, []);
axis image;
title('"Keeper" blobs');

% Now use the keeper blobs as a mask on the original image.
% This will let us display the original image in the regions of the keeper blobs.
maskedImage = img2; % Simply a copy at first.
maskedImage(~keeperBlobsImage) = 0;  % Set all non-keeper pixels to zero.
subplot(3, 3, 8);
imshow(maskedImage);
axis image;
title('');



%Begin Analyzation of Right Half

ImgRight = imread('rightHalf_Brain.jpeg');
figure
imshow(ImgRight);

%Converts image to Black and white

img2=im2bw(ImgRight,graythresh(ImgRight));
imshow(img2)
%Changes black to white and vice versa in image
img2=~img2;
imshow(img2)

%Finds the boundaries for each object and stores in variable B
%Number of Objects found is printed on Image by Text Function  
B = bwboundaries(img2);
imshow(img2)
text(4,4,strcat('\color{green}Objects Found:',num2str(length(B))))
hold on 

%Marks the boundaries of objects with green lines
for cnt = 1:length(B)
boundary = B{cnt};
plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1)
end

imshow(img2)

% Identify individual blobs by seeing which pixels are connected to each other.
labelednImage = bwlabel(img2, 8);     % Label each blob so we can make measurements of it
% labeledImage is an integer-valued image where all pixels in the blobs have values of 1, or 2, or 3, or ... etc.
subplot(3, 3, 1);
imshow(labelednImage, []);  % Show the gray scale image.

title('Labeled Image, from bwlabel()');

% Let's assign each blob a different color to visually show the user the distinct blobs.
coloredLabels = label2rgb (labelednImage, 'hsv', 'k', 'shuffle'); % pseudo random color labels
% coloredLabels is an RGB image.  We could have applied a colormap instead
subplot(3, 3, 5);
imshow(coloredLabels);

axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
caption = sprintf('Pseudo colored labels, from label2rgb().\nBlobs are numbered from top to bottom, then from left to right.');
title(caption);

% Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
blobMeasurements = regionprops(labelednImage, img2, 'all');
numberOfBlobs = size(blobMeasurements, 1);

% bwboundaries() returns a cell array, where each cell contains the row/column coordinates for an object in the image.
% Plot the borders of all the coins on the original grayscale image using the coordinates returned by bwboundaries.
subplot(3, 3, 6);
imshow(img2);
title('Outlines, from bwboundaries()'); 
axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
hold on;
boundaries = bwboundaries(img2,4);
numberOfBoundaries = size(boundaries, 1);
for cnt = 1 : numberOfBoundaries
	thisBoundary = boundaries{cnt};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
end
hold off;

textFontSize = 14;	% Used to control size of "blob number" labels put atop the image.
labelShiftX = -7;	% Used to align the labels in the centers of the coins.
blobECD = zeros(1, numberOfBlobs);

% Loop over all blobs printing their measurements to the command window.
for count = 1 : numberOfBlobs           % Loop through all blobs.
	% Find the mean of each blob.  
	% directly into regionprops.  The way below works for all versions including earlier versions.)
	thisBlobsPixels = blobMeasurements(count).PixelIdxList;  % Get list of pixels in current blob.
	meanGL = mean(img2(thisBlobsPixels)); % Find mean intensity (in original image!)
	meanGL2008a = blobMeasurements(count).MeanIntensity; % Mean again, but only for version >= R2008a
	
	blobArea = blobMeasurements(count).Area;		% Get area.
	blobPerimeter = blobMeasurements(count).Perimeter;		% Get perimeter.
	blobCentroid = blobMeasurements(count).Centroid;		% Get centroid one at a time
	blobECD(count) = sqrt(4 * blobArea / pi);					% Compute ECD - Equivalent Circular Diameter.
    if blobArea > 5000
    fprintf('The Area of the Right Half: %d\n', blobArea);
    healthymeasurements(2,:) = [2 meanGL blobArea blobPerimeter blobECD(count)]; %Matrix to store collected data
    end
    
    % Put the "blob number" labels on the "boundaries" grayscale image.
    text(blobCentroid(1) + labelShiftX, blobCentroid(2), num2str(count));
end
xlswrite('Export_GajavelliData.xlsx',healthymeasurements, 'Data', 'A2')



allBlobAreas = [blobMeasurements.Area];
% Get a list of the blobs that meet our criteria and we need to keep.
allowableAreaIndexes = allBlobAreas > 2000; % Take the small objects.
keeperIndexes = find(allowableAreaIndexes);
% Extract only those blobs that meet our criteria, and
% eliminate those blobs that don't meet our criteria.
% Note how we use ismember() to do this.  Result will be an image - the same as labeledImage but with only the blobs listed in keeperIndexes in it.
keeperBlobsImage = ismember(labelednImage, keeperIndexes);
% Re-label with only the keeper blobs kept.
labelednImage = bwlabel(keeperBlobsImage, 8);     % Label each blob so we can make measurements of it
subplot(3, 3, 7);
imshow(labelednImage, []);
axis image;
title('"Keeper" blobs');

% Now use the keeper blobs as a mask on the original image.
% This will let us display the original image in the regions of the keeper blobs.
maskedImage = img2; % Simply a copy at first.
maskedImage(~keeperBlobsImage) = 0;  % Set all non-keeper pixels to zero.
subplot(3, 3, 8);
imshow(maskedImage);
axis image;
title('');
  

    figure 
    montage({I,ImgChart,imgGray,ImgLeft,ImgRight})
    multi =  {I,ImgChart,imgGray,ImgLeft,ImgRight};
    z = montage(multi,'Size',[2,3]);
    montage_IM = z.CData;
    %write to file
    
    stringName = strcat(fullFileName, '_Montage.jpeg');
    baseFileName = sprintf(stringName);
    fullFileName2 = fullfile('Montages', baseFileName);
    imwrite(montage_IM,fullFileName2);

  
         end
end


