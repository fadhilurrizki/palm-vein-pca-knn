%===============awal==========
% image = imread('dataset/training/001_r_940_01.jpg');
% image = uint8(image);
% image = imrotate(image,270);
% %image = imcrop(image,[10 40 500 580])
% image = imnoise(image,'salt & pepper',0.02);
% image = medfilt2(image,[3 3]);
% level = graythresh(image);
% image = im2bw(image,level);
% image = bwareaopen(image,3000);
% a = regionprops(image,'Image');
% roi = a.Image;
% figure,imshow(roi);
%================================
orang = '';
jumlahorang = 3;
jumlahcitra = 3;
D = zeros([10000 jumlahcitra]);
counter = 0;
    for i = 1:3
        filename =  strcat('training/001_r_940_0',int2str(i),'.jpg');
        img = imread(filename);
        img = imrotate(img,270);
        im2=img;
        img = imnoise(img,'salt & pepper',0.02);
        img = medfilt2(img,[3 3]);
        level = graythresh(img);
        img = im2bw(img,level);
        img = bwareaopen(img,3000);
        a = regionprops(img,'Image','BoundingBox');
        Box = a.BoundingBox;
        img = imcrop(img,Box);
%         img = imresize(img,[100 100]);
%         [maxtab, mintab] = peakdet(img, 0.5);
%         hold on; plot(mintab(:,1), mintab(:,2), 'g*');
%         figure, plot(maxtab(:,1), maxtab(:,2), 'r*');
        x = 1;
        counter = counter +1;
        for k=1:100
            for j=1:100
                D(x,counter) = img(k,j);
                x=x+1;
            end
        end
    end

% D = {d1,d2,...,dn}, n = 10000 (dimensi citra)
meann = zeros([10000 1]);

    for j=1:10000
        meann(j,1) = mean(D(j,:));
    end

%D - mean
Y = zeros([10000 jumlahcitra]);
for i=1:3
    Y(:,i) = D(:,i)-meann;
end
%Covariance
A = transpose(Y)*Y;
%nilai eigen & vektor eigen
[v,d] = eig(A);
[d order] = sort(diag(d), 'descend');
v = v(:,order);
pctrain = Y*v;

jumlahcitra = 2;
D = zeros([10000 jumlahcitra]);
counter = 0;
    for i = 1:jumlahcitra
        filename =  strcat('testing/001_r_940_0',int2str(i+4),'.jpg');
        img = imread(filename);
        img = imrotate(img,270);
        img = imnoise(img,'salt & pepper',0.02);
        img = medfilt2(img,[3 3]);
        level = graythresh(img);
        img = im2bw(img,level);
        img = bwareaopen(img,3000);
        a = regionprops(img,'Image','BoundingBox');
        Box = a.BoundingBox;
        img = imcrop(img,Box);
        img = imresize(img,[100 100]);
        x = 1;
        counter = counter +1;
        for k=1:100
            for j=1:100
                D(x,counter) = img(k,j);
                x=x+1;
            end
        end
    end

% D = {d1,d2,...,dn}, n = 10000 (dimensi citra)
meann = zeros([10000 1]);

    for j=1:10000
        meann(j,1) = mean(D(j,:));
    end

%D - mean
Y = zeros([10000 jumlahcitra]);
for i=1:jumlahcitra
    Y(:,i) = D(:,i)-meann;
end
%Covariance
A = transpose(Y)*Y;
%nilai eigen & vektor eigen
[v,d] = eig(A);
[d order] = sort(diag(d), 'descend');
v = v(:,order);
pctes = Y*v;
% md1 = ClassificationKNN.fit(pc,pc2);
% disp(pc2);
% W = transpose(pc)*Y;
% disp(W);
pcaa = pca(D);
[coeff score latent] = pca(D,'Algorithm','eig');
disp(coeff);
% disp(score);
% disp(latent);

Sample = coeff;
Training=[0 0;.654 .343;.3 .4];
Group = [1;2;3];
disp(Training);
Class = knnclassify(Sample, Training, Group)

% A = zeros([10000 2]);
% for i = 1:2
%     filename =  strcat('dataset/testing/002_r_940_0',int2str(i+4),'.jpg');
%     img = imread(filename);
%     img = imrotate(img,270);
%     img = imnoise(img,'salt & pepper',0.02);
%     img = medfilt2(img,[3 3]);
%     level = graythresh(img);
%     img = im2bw(img,level);
%     img = bwareaopen(img,3000);
%     a = regionprops(img,'Image','BoundingBox');
%     Box = a.BoundingBox;
%     img = imcrop(img,Box);
%     img = imresize(img,[100 100]);
%     x = 1;
%     for k=1:100
%         for j=1:100
%             A(x,i) = img(k,j);
%             x=x+1;
%         end
%     end
% end
% 
% pc = pca(A);
% disp(pc);