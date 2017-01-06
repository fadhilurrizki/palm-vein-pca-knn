function varargout = palmveinmain(varargin)

% PALMVEINMAIN MATLAB code for palmveinmain.fig
%      PALMVEINMAIN, by itself, creates a new PALMVEINMAIN or raises the existing
%      singleton*.
%
%      H = PALMVEINMAIN returns the handle to a new PALMVEINMAIN or the handle to
%      the existing singleton*.
%
%      PALMVEINMAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PALMVEINMAIN.M with the given input arguments.
%
%      PALMVEINMAIN('Property','Value',...) creates a new PALMVEINMAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before palmveinmain_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to palmveinmain_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help palmveinmain

% Last Modified by GUIDE v2.5 07-May-2015 07:54:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @palmveinmain_OpeningFcn, ...
                   'gui_OutputFcn',  @palmveinmain_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before palmveinmain is made visible.
function palmveinmain_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to palmveinmain (see VARARGIN)

% Choose default command line output for palmveinmain
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes palmveinmain wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = palmveinmain_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in training.
function training_Callback(hObject, eventdata, handles)
% hObject    handle to training (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in testing.
function testing_Callback(hObject, eventdata, handles)
% hObject    handle to testing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global IMG;
global data_latih;
global kelas_latih;
T = zeros([6400 1]);
%         filename =  strcat('training/020_r_940_01.jpg');
%         img = imread(filename);
         img = IMG;
        img = imrotate(img,270);
        img = imnoise(img,'salt & pepper',0.02);
        img = medfilt2(img,[3 3]);   
        im2=img;
        level = graythresh(img);
        img = im2bw(img,level);
        img = bwareaopen(img,50000);
        a = regionprops(img,'Centroid')
        centx = ceil(a.Centroid(1));
        centy = ceil(a.Centroid(2));
        img = imcrop(im2,[centx-40 centy-40 80 80]);
        img = histeq(img);
        x = 0;
        counter = 1;
        for k=1:80
            for j=1:80
                x=x+1;
                T(x,counter) = img(k,j);
            end
        end
        axes(handles.hasil);
        imshow(img);
k = 3;
kelas_uji = knnclassify(T', data_latih', kelas_latih', k);
s = kelas_uji;
set(handles.hasil2,'string',s);



% global pc;
% global vec;
% t = zeros([6400 1]);
% filename =  'testing/002_r_940_05.jpg';
%         img = imread(filename);
%         img = imrotate(img,270);
%         img = imnoise(img,'salt & pepper',0.02);
%         img = medfilt2(img,[3 3]);
%         
%          im2=img;
%         level = graythresh(img);
%         img = im2bw(img,level);
%         img = bwareaopen(img,50000);
%         a = regionprops(img,'Centroid')
%         centx = ceil(a.Centroid(1));
%         centy = ceil(a.Centroid(2));
%         img = imcrop(im2,[centx-40 centy-40 80 80]);
%         img = histeq(img);
%         x = 1;
%         counter =1;
%         for k=1:80
%             for j=1:80
%                 t(x,counter) = img(k,j);
%                 x=x+1;
%             end
%         end
%         axis(handles.hasil)
%         imshow(img);
% %t=t-m;
% %Project the testing data on the space of the training data
% t=vec'*t;
% 
% % Then if you want to know what is the class of this test data?  just use
% % any classifier (In our case we used minimum distance classifier)
% alldata=t';
% alldata(2:size(pc,2)+1,:)=pc';
% dist=pdist(alldata);
% [a,b]=min(dist(:,1:size(x,2)));
% disp(b);



% % formku=guidata(gcbo);
% % I=get(formku.namaFile, 'String');
% % ROI1 = findROI(I); % Call function
% % ROI1 = histeq(ROI1);
% % sizeroi = size(ROI1);
% % x = zeros([sizeroi(1)*sizeroi(2) 2]);
% % a=0;
% % for i=1:sizeroi(1)
% %     for j=1:sizeroi(2)
% %         a=a+1;
% %         x(a,1) = ROI1(i,j);
% %         x(a,2) = ROI2(i,j);
% %     end
% % end
% % coef = pca(x);
% formku=guidata(gcbo);
% jumlahcitra = 2;
% D = zeros([10000 jumlahcitra]);
% counter = 0;
%     for i = 1:jumlahcitra
%         filename =  strcat('testing/002_r_940_0',int2str(i+4),'.jpg');
%         img = imread(filename);
%         img = imrotate(img,270);
%         img = imnoise(img,'salt & pepper',0.02);
%         img = medfilt2(img,[3 3]);
%         im2=img;
%         level = graythresh(img);
%         img = im2bw(img,level);
%         img = bwareaopen(img,3000);
% %         a = regionprops(img,'Image','BoundingBox');
% %         Box = a.BoundingBox;
% %         img = imcrop(img,Box);
%         a = regionprops(img, 'Image', 'Centroid');
%         centx = ceil(a.Centroid(1));
%         centy = ceil(a.Centroid(2));
%         img = imcrop(im2,[centx-50 centy-50 100 100]);
%         imshow(img);
% %          ROI(img == 1) = 1;
% %          ROI(img ~= 1) = 0;
% %          SEG = image.*ROI;
% %          img = SEG;
% 
%         img = imresize(img,[100 100]);
%         axes(handles.hasil)
%         imshow(img);
%         x = 1;
%         counter = counter +1;
%         for k=1:100
%             for j=1:100
%                 D(x,counter) = img(k,j);
%                 x=x+1;
%             end
%         end
%     end
% 
% % D = {d1,d2,...,dn}, n = 10000 (dimensi citra)
% meann = mean(D')';
% 
% %D - mean
% Y = D - repmat(meann,1,jumlahcitra);
% 
% %Covariance
% A = transpose(Y)*Y;
% %nilai eigen & vektor eigen
% [v,d] = eig(A);
% [d order] = sort(diag(d), 'descend');
% v = v(:,order);
% pctes = Y*v;
% % md1 = ClassificationKNN.fit(pc,pc2);
% % disp(pc2);
% % W = transpose(pc)*Y;
% % disp(W);
% pcaa = pca(D);
% [coeff score latent] = pca(D,'Algorithm','eig');
% % disp(coeff);
% % disp(score);
% % disp(latent);
% 
% Sample = coeff;
% Training=[0 0;.654 .343;.3 .4];
% Group = [1;2;3];
% % disp(Training);
% Class = knnclassify(Sample, Training, Group);
% 
% 
% h1=strcat(Class);
% set(formku.hasilklasifikasi,'String',num2str(Class));
% disp(Class);

% If you have test data do the following
%t=[1;1]  % this test data is close to the first class
%Subtract the mean from the test data


function hasilklasifikasi_Callback(hObject, eventdata, handles)
% hObject    handle to hasilklasifikasi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hasilklasifikasi as text
%        str2double(get(hObject,'String')) returns contents of hasilklasifikasi as a double


% --- Executes during object creation, after setting all properties.
function hasilklasifikasi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hasilklasifikasi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in open.
function open_Callback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global IMG;
[namafile, formatfile] = uigetfile({'*.jpg'}, 'membuka gambar'); %memilih gambar
IMG = imread([formatfile, namafile]); %membaca gambar
IMG=uint8(IMG);
handles.filename = namafile;
guidata(hObject, handles);
axes(handles.hasil); %memilih axes1 sebagai letak gambar yang dimunculkan
imshow(IMG);


function namaFile_Callback(hObject, eventdata, handles)
% hObject    handle to namaFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of namaFile as text
%        str2double(get(hObject,'String')) returns contents of namaFile as a double


% --- Executes during object creation, after setting all properties.
function namaFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to namaFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hasil2_Callback(hObject, eventdata, handles)
% hObject    handle to hasil2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hasil2 as text
%        str2double(get(hObject,'String')) returns contents of hasil2 as a double


% --- Executes during object creation, after setting all properties.
function hasil2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hasil2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadtraining.
function loadtraining_Callback(hObject, eventdata, handles)
% hObject    handle to loadtraining (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global vec2;
global pc;
global data_latih;
global kelas_latih;
global projectedtraining;
jumlahorang = 100;
orang = '';
jumlahcitra = 4;
D = zeros([6401 jumlahorang*jumlahcitra]);
counter = 0;
jumlahdata = jumlahcitra*jumlahorang;
nocitra = '';
nocitra = [1 2 5 6];
for l=1:jumlahorang
    for i = 1:jumlahcitra
        if(l<10)
            orang = strcat('00',int2str(l));
        else
        if(l<100)
            orang = strcat('0',int2str(l));
        else
            orang = int2str(l);
        end
        end
        
        filename =  strcat('training/',orang,'_r_940_0',int2str(nocitra(i)),'.jpg');
        img = imread(filename);
        img = imrotate(img,270);
        img = imnoise(img,'salt & pepper',0.02);
        img = medfilt2(img,[3 3]);
        
         im2=img;
        level = graythresh(img);
        img = im2bw(img,level);
        img = bwareaopen(img,50000);
        a = regionprops(img,'Centroid')
        centx = ceil(a.Centroid(1));
        centy = ceil(a.Centroid(2));
        img = imcrop(im2,[centx-40 centy-40 80 80]);
        img = histeq(img);
        x = 1;
        counter = counter +1;
        D(1,counter) = l;
        for k=1:80
            for j=1:80
                x=x+1;
                D(x,counter) = img(k,j);
            end
        end
        
    end
end
data_latih = D(2:6401,:);
kelas_latih = D(1,:);

s = 'Finish';
set(handles.teks,'string',s);

% [r,c]=size(data_latih);
% % Compute the mean of the data matrix "The mean of each row"
% m=mean(data_latih')';
% % Subtract the mean from each image [Centering the data]
% d=data_latih-repmat(m,1,c);
% 
% % Compute the covariance matrix (co)
% co=d*d';
% 
% %Compute the eigen values and eigen vectors of the covariance matrix
% [eigvector,eigvl]=eig(co);
% 
% % Sort the eigen vectors according to the eigen values
% eigvalue = diag(eigvl);
% [junk, index] = sort(-eigvalue);
% eigvalue = eigvalue(index);
% eigvector = eigvector(:, index);
% 
% count1=0;
% for i=1:size(eigvalue,1)
%     if(eigvalue(i)>0)
%         count1=count1+1;
%     end
% end
% 
% vec=eigvector(:,1:count1);
% 
% pc=vec'*d;
% disp(pc);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_latih;
global kelas_latih;
global vec2;
global pc;
kelasuji = zeros([200 2]);
jumlahorang = 100;
jumlahcitra = 2;
counter = 0;
nocitra = [3 6];
T = zeros([6400 200]);
for l=1:jumlahorang
    for i = 1:jumlahcitra
        if(l<10)
            orang = strcat('00',int2str(l));
        else
        if(l<100)
            orang = strcat('0',int2str(l));
        else
            orang = int2str(l);
        end
        end
        filename =  strcat('training/',orang,'_r_940_0',int2str(nocitra(i)),'.jpg');
        img = imread(filename);
        img = imrotate(img,270);
        img = imnoise(img,'salt & pepper',0.02);
        img = medfilt2(img,[3 3]);   
        im2=img;
        level = graythresh(img);
        img = im2bw(img,level);
        img = bwareaopen(img,50000);
        a = regionprops(img,'Centroid')
        centx = ceil(a.Centroid(1));
        centy = ceil(a.Centroid(2));
        img = imcrop(im2,[centx-40 centy-40 80 80]);
        img = histeq(img);
        x = 0;
        counter = counter + 1;
        for k=1:80
            for j=1:80
                x=x+1;
                T(x,counter) = img(k,j);
            end
        end     
        kelasuji(counter,1) = l;
    end
end
[r,c]=size(T);
% Compute the mean of the data matrix "The mean of each row"
m=mean(T')';
d=T-repmat(m,1,c);
pc2=vec2'*d;
k = 1;
for i=1:200
    kelas_uji = knnclassify(pc2(:,i)', pc', kelas_latih', k);
    kelasuji(i,2) = kelas_uji;
    set(handles.uitable1,'Data',kelasuji);
end
jlh = 0;
for i=1:200
    if(kelasuji(i,1)==kelasuji(i,2))
        jlh = jlh+1;
    end
end
akurasi = jlh/200*100;
s = strcat('Akurasi = ',num2str(akurasi),' %');
set(handles.hasilklasifikasi,'string',s);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_latih;
global kelas_latih;
global pc;
global vec2;
[r,c]=size(data_latih);
% Compute the mean of the data matrix "The mean of each row"
m=mean(data_latih')';
% Subtract the mean from each image [Centering the data]
d=data_latih-repmat(m,1,c);

% Compute the covariance matrix (co)
co=d'*d;

%Compute the eigen values and eigen vectors of the covariance matrix
[eigvector,eigvl]=eig(co);

% Sort the eigen vectors according to the eigen values
eigvalue = diag(eigvl);
[junk, index] = sort(-eigvalue);
eigvalue = eigvalue(index);
eigvector = eigvector(:, index);

count1=0;
for i=1:size(eigvalue,1)
    if(eigvalue(i)>0)
        count1=count1+1;
    end
end

vec=eigvector(:,1:150);
vec2 = d*vec;

pc=vec2'*d;
size(pc);
s='Finish';
set(handles.text3,'string',s),