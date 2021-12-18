% READ FILE
clc; clear all; close all;
test_image = imread('40_training.tif');
ref_image = imread('40_training.png');
green_channel_image = test_image(:,:,2);
background_image = green_channel_image;
background_image = im2bw(background_image, 20/255);
background_image = ~background_image;

% CLAHE
enhanced_image = adapthisteq(green_channel_image,'NumTiles',[75 75],'ClipLimit',0.09);

% FILTER
filtered_image = imgaussfilt(enhanced_image, 0.1);
filtered_image = medfilt2(filtered_image);
filtered_image = imadjust(filtered_image, [60/255 127/255], [0.1 1]);

% ISODATA TRHESHOLDING
threshold_image = im2double(filtered_image);
level = threshold_level(threshold_image);
threshold_image = imbinarize(threshold_image, level);

% MORPHOLOGICAL
vessel_image = ~threshold_image;
vessel_image = imsubtract(vessel_image, background_image);
figure, imshow(vessel_image)
vessel_image = bwareaopen(vessel_image, 150);
vessel_image = imsubtract(vessel_image, background_image);
vessel_image = medfilt2(vessel_image);
[row, col] = size(vessel_image);
for i = 1:row
    for j = 1:col
        if(vessel_image(i,j)<0)
            vessel_image(i,j) = 0;
        end
    end
end
figure, imshow(vessel_image)
% title('Segmentasi Pembuluh Darah Retina')

% PERFOEMANCE ALGORITHM
[TP, FP, TN, FN] = performance(vessel_image, ref_image) 

vessel_image = imfuse(vessel_image, ref_image, 'falsecolor');
figure, imshow(vessel_image);
% LEVEL ISODATA THRESHOLDING
function level = threshold_level(image)
    image = im2uint8(image(:));
    [histogram_count, bin_number] = imhist(image);
    i = 1;
    cumulative_sum = cumsum(histogram_count);
    T(i) = (sum(bin_number.*histogram_count))/cumulative_sum(end);
    T(i) = round(T(i));
    
    cumulative_sum_2 = cumsum(histogram_count(1:T(i)));
    MBT = sum(bin_number(1:T(i)).*histogram_count(1:T(i)))/cumulative_sum_2(end);
    
    cumulative_sum_3 = cumsum(histogram_count(T(i):end));
    MAT = sum(bin_number(T(i):end).*histogram_count(T(i):end))/cumulative_sum_3(end);
    
    i = i+1;
    T(i) = round((MAT+MBT)/2);
    
    while abs(T(i)-T(i-1))>=1
        cumulative_sum_2 = cumsum(histogram_count(1:T(i)));
        MBT = sum(bin_number(1:T(i)).*histogram_count(1:T(i)))/cumulative_sum_2(end);

        cumulative_sum_3 = cumsum(histogram_count(T(i):end));
        MAT = sum(bin_number(T(i):end).*histogram_count(T(i):end))/cumulative_sum_3(end);

        i = i+1;
        T(i) = round((MAT+MBT)/2);
        
        threshold = T(i);
    end
    
    level = (threshold-1)/(bin_number(end)-1);
end
    
function [TP, FP, TN, FN] = performance(algo_image, ref_image) 
    [row, col] = size(algo_image);
    TP=0; FP=0; TN=0; FN=0;
    for i = 1:row
        for j = 1:col
            if((algo_image(i,j)==1) && (ref_image(i,j)==255));
                TP = TP+1;
            elseif((algo_image(i,j)==1) && (ref_image(i,j)==0));
                FP = FP+1;
            elseif((algo_image(i,j)==0) && (ref_image(i,j)==255));
                FN = FN+1;
            else((algo_image(i,j)==0) && (ref_image(i,j)==0));
                TN = TN+1;
            end
        end
    end
end
