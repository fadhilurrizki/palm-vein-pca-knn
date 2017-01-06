
orang = '';
jumlahorang = 3;
jumlahcitra = 3;
counter = 0;
D = zeros([10000 jumlahcitra]);
img = imread('training/001_r_940_01.jpg');
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
img = imresize(img,[100 100]);
[maxtab, mintab] = peakdet(img, 0.1);

x = 1;
counter = counter +1;
for k=1:100
    for j=1:100
        D(x,counter) = img(k,j);
        x=x+1;
    end
end
meann = zeros([10000 1]);
meann(j,1) = mean(D(j,:));
Y = zeros([10000 jumlahcitra]);
Y(:,1) = D(:,1)-meann;
A = transpose(Y)*Y;
%nilai eigen & vektor eigen
[v,d] = eig(A);
[d order] = sort(diag(d), 'descend');
v = v(:,order);
pc = Y*v;
figure, imshow(img);
