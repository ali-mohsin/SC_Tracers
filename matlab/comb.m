clc;
clear;
close all;

load('image_pm1.mat');

image_d = imread('transparent_soft_shadow.ppm','ppm');

image_d = double(image_d)/255;

image = image_d+image_c+image_g;

imshow(uint8(image*255));