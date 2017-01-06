function ROIout = findROI(imname)
A = imread(imname);
A=imrotate(A,270);
% A = imname;
G = fspecial('gaussian',[5 5], 2);      % Applying the Smoothing Filter
C = imfilter(A, G, 'same');
[w, h] = size(C);


BW = im2bw(C, graythresh(C));           % Convert output into Binary
[B,~,~] = bwboundaries(BW);             % Detect the Boundary

centroid = regionprops(BW,'Centroid');  % Locate the centroid
centroid = centroid.Centroid;
% Cendroid(1) - X coordinate of Centroid
% Cendroid(2) - Y coordinate of Centroid

outline = flipud(B{1});     % This is the edge we are interested in

regionLine = zeros(length(outline),3);  % Put all coordinates and distances in the same array

for i = 1:length(outline)
    regionLine(i,:) = [outline(i,2) outline(i,1) sqrt((outline(i,2)-centroid(1))^2+(outline(i,1)-centroid(2))^2)];
end

% Use -regionLine (flip the line up side down) to find minima...
[~, i] = findpeaks(smooth(-regionLine(:,3),50),'MINPEAKDISTANCE',10,'MINPEAKHEIGHT',mean(-regionLine(:,3)));
% Higher order Moving average filter used... User may wish to change the
% parameters of the above function

if numel(regexp(imname,'_l_')) %Left Image minima of interest 1:2 and end-1:end
    coordinates = [regionLine(i(1:2),1:2); regionLine(i(end-1:end),1:2)];
else %Right Image minima of interest 1:3 and end
    coordinates = [regionLine(i(1:3),1:2); regionLine(i(end),1:2)];
end

sortCoord = sortrows(coordinates,2);    % Sort the coordinates from top of image to bottom

sc = lineScore(sortCoord);              % Distance based score to reject one point
if sc(1)<sc(2)
    sortCoord(4,:) = [];                % Reject point 4 (1, 2, 3 are close)
else
    sortCoord(1,:) = [];                % Reject point 1 (2, 3, 4 are close)
end

% Keeping track of the centroid after rotation
zeroIm1 = zeros(w,h);
zeroIm2 = zeros(w,h);

zeroIm1(round(sortCoord(1,2))-2:round(sortCoord(1,2))+2,round(sortCoord(1,1))-2:round(sortCoord(1,1))+2) = 255;
zeroIm2(round(sortCoord(3,2))-2:round(sortCoord(3,2))+2,round(sortCoord(3,1))-2:round(sortCoord(3,1))+2) = 255;

% Angle of rotation to keep point 1 and 3 in the vertical...
grad = (atan((sortCoord(3,2)-sortCoord(1,2))/(sortCoord(3,1)-sortCoord(1,1))))*180/pi;

if numel(regexp(imname,'_l_')) %Left Image
    rotPalm = imrotate(A,grad+90);          % Rotation Angle
    rotzeroIm1 = imrotate(zeroIm1,grad+90);
    rotzeroIm2 = imrotate(zeroIm2,grad+90);
else
    rotPalm = imrotate(A,grad-90);          % Rotation Angle
    rotzeroIm1 = imrotate(zeroIm1,grad-90);
    rotzeroIm2 = imrotate(zeroIm2,grad-90);
end

[top2, ~] = find(rotzeroIm1==255);          % Find Top coordinate in Rotated Image
[bottom2, ~] = find(rotzeroIm2==255);       % Find Bottom coordinate in Rotated Image


xOffset = round(centroid(1));      % This is the xaxis value along the bisector...
 % Should be tuned...

new_c = [xOffset round((top2(1)+bottom2(1))/2)];

ROIout = rotPalm(new_c(2)-63:new_c(2)+64,new_c(1)-63:new_c(1)+64);

%% Showing an output - Comment/Uncomment as needed
%figure
% subplot 121
% imshow(rotPalm)
% hold on
% plot(new_c(1),new_c(2),'b*')
% rectangle('Position',[new_c(1)-63,new_c(2)-63,128,128])

% title('Input Image')
% 
% subplot 122
% imshow(ROIout)
% title('Region of Interest')
