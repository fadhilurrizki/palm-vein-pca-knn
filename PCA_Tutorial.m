% prophet mohammed said [ALLAH will help any one helped his/her brother/sister] PBUH

%This code to apply PCA (Principal Component Analysis) 
% for any information please send to engalaatharwat@hotmail.com
%Egypt - HICIT - +20106091638


% Remember that each column of the data matrix(input matrix) represent one image or pattern  
% Note: the data here represent two classes
% Class 1: data(:,1:4)
% Class 2: data(:,5:8)

% data = [1     1     2     0     7     6     7     8
%         3     2     3     3     4     5     5     4];
jumlahorang = 4;
jumlahcitra = 4;
D = zeros([400 jumlahorang*jumlahcitra]);
counter = 0;
jumlahdata = jumlahcitra*jumlahorang;
for l=1:jumlahorang
    for i = 1:jumlahcitra
        filename =  strcat('training/00',int2str(l),'_r_940_0',int2str(i),'.jpg');
        img = imread(filename);
        img = imrotate(img,270);
        img = imnoise(img,'salt & pepper',0.02);
        img = medfilt2(img,[3 3]);
        
         im2=img;
        level = graythresh(img);
        img = im2bw(img,level);
        img = bwareaopen(img,10000);
        a = regionprops(img,'Centroid')
        centx = ceil(a.Centroid(1));
        centy = ceil(a.Centroid(2));
        img = imcrop(im2,[centx-50 centy-50 100 100]);
        img = imresize(img,[20 20]);
        x = 1;
        counter = counter +1;
        for k=1:20
            for j=1:20
                D(x,counter) = img(k,j);
                x=x+1;
            end
        end
    end
end

[r,c]=size(D);
% Compute the mean of the data matrix "The mean of each row"
m=mean(D')';
% Subtract the mean from each image [Centering the data]
d=D-repmat(m,1,c);



% Compute the covariance matrix (co)
co=d*d';

%Compute the eigen values and eigen vectors of the covariance matrix
[eigvector,eigvl]=eig(co);

% Sort the eigen vectors according to the eigen values
eigvalue = diag(eigvl);
[junk, index] = sort(-eigvalue);
eigvalue = eigvalue(index);
eigvector = eigvector(:, index);

% Compute the number of eigen values that greater than zero (you can select any threshold)
count1=0;
for i=1:size(eigvalue,1)
    if(eigvalue(i)>0)
        count1=count1+1;
    end
end

vec=eigvector(:,1:count1);

% Compute the feature matrix (the space that will use it to project the testing image on it)
x=vec'*d;
disp(x);


% If you have test data do the following
%t=[1;1]  % this test data is close to the first class
%Subtract the mean from the test data
t = zeros([400 1]);
filename =  'testing/002_r_940_05.jpg';
        img = imread(filename);
        img = imrotate(img,270);
        img = imnoise(img,'salt & pepper',0.02);
        img = medfilt2(img,[3 3]);
        
         im2=img;
        level = graythresh(img);
        img = im2bw(img,level);
        img = bwareaopen(img,10000);
        a = regionprops(img,'Centroid');
        figure,imshow(img);
        centx = ceil(a.Centroid(1));
        centy = ceil(a.Centroid(2));
        img = imcrop(im2,[centx-10 centy-10 20 20]);

        x = 1;
        counter =1;
        for k=1:20
            for j=1:20
                t(x,counter) = img(k,j);
                x=x+1;
            end
        end

t=t-m;
%Project the testing data on the space of the training data
t=vec'*t;

% Then if you want to know what is the class of this test data?  just use
% any classifier (In our case we used minimum distance classifier)
alldata=t';
disp(alldata);
alldata(2:size(x,2)+1,:)=x';
dist=pdist(alldata);
[a,b]=min(dist(:,1:size(x,2)));
disp(b);
% So, b determine the closest observation to test data;