clc
close all;
clear all;
%[filename, pathname] = uigetfile({'*.*';'*.bmp';'*.jpg';'*.png'}, 'Pick a Image File');
%myImage = imread([pathname,filename]);


fileFolder = 'F:\BTP_2\project\test'; 
filePattern = fullfile(fileFolder, '*.png');
dirOutput = dir(filePattern);

fileNames = {dirOutput.name};
numberOfImageFiles = numel(fileNames);
str = {};
load('dataset.mat');
for j= 1:numberOfImageFiles
    myImage = imread(fileNames{j});
    
    J = rgb2gray(myImage);
    img = im2bw(J, 0.9);
        
   
    
    [H, W] = size(img);
	A = W*H;
	T = 0;
     % calculate the total number of black pixels in the images
    for x =1:W
		for y =1:H
			if img(y,x) == 0
				T = T +1;
            end
        end
        
    end	
      % store the width,height, area and T(number balck pixels) in a list
 basicGlobalFeature = [W H A T];
 
 
	Ci = (4 * W * H) / (3.14 * (W*W + H*H));
	Srad = sqrt(W*W + H*H)/2.0;
 circularityFeature = [Ci,Srad];
        
 
 % calculate the number of black pixels in each row and store in a list
     Pv = [];
    for y =1:H
		tB = 0;
		for x=1:W
			if img(y, x) == 0 
				tB = tB +1;
            end
        end
		Pv(y)= tB;
    end
    
    
   % calculate the number of black pixels in each column and store in a list
   Ph = [];
    for x =1:W
		tB = 0;
		for y=1:H
			if img(y, x) == 0 
				tB = tB +1;
            end
        end
		Ph(x)= tB;
    end
	
    % finding a vertical center. it is just finding center of mass of
    % object
    total = 0;
	for y =1:H
		total = total + y * Pv(y);
    end
	verticalCenter = total/T;
   
    % horzontal center
    total = 0;
    for x =1:W
		total = total + x* Ph(x);
    end
	horizontalCenter = total/T;
    
    
    % the_y is a row index having most number of black pixels among all the rows and the_value is the value at that index
    the_y = 0;
	the_value = 0;
    for y =1:H
		if Pv(y) > the_value
			the_value = Pv(y);
			the_y = y;
        end
    end
	globalbaseLine = [the_y, the_value];
    
    
    BL = the_y;
    value = the_value;
	upperLimit = 0;
	diff = 0;
	for y =1:BL
		smoothV = y * value / BL;
		tempDiff = abs(Pv(y) - smoothV);
		if tempDiff > diff
			upperLimit = y;
			diff = tempDiff;
        end
    end
    
    lowerLimit = 0;
    diff = 0;
    for y =BL:H
		smoothV = value - (y * value / ((H-1) - BL));
		tempDiff = abs(Pv(y) - smoothV);
		if tempDiff > diff
			lowerLimit = y;
			diff = tempDiff;
        end
    end
    
    
   % creating some features from the above properties 
    HtW = (H) / W;
	AtC = Ci;
	TtA = (T) / A;
	BtH = (the_y) / H;
	LtH = (lowerLimit) / H;
	UtH = (H-upperLimit+1) / H;
    
    
    testarray = zeros(1,8);
    i=1;
    
    % storing the features in the list
    testarray(i,1) =HtW;
    testarray(i,2) =AtC;
   
    testarray(i,3) = TtA;
    testarray(i,4) =BtH;

    testarray(i,5) = LtH;
    testarray(i,6) = UtH;
    testarray(i,7) = verticalCenter;
    testarray(i,8) = horizontalCenter;
  
    
    name = {'jack';'elliot';'tom';'swap'};


result = knnclassify(testarray,centroid,name);
 
 str{j} = result;
 %disp(result);

end

for k=1:numberOfImageFiles
    disp(str{k});
    
end



