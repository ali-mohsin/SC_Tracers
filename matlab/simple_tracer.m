clc;
clear;
close all;

width = 128;
height = 128;

image = zeros(height,width,3);

FOV = 40;
aspectratio = width/height;

angle = tand(0.5*FOV);

fin = fopen('prism31.asc','r');
%points = fscanf(fin,'triangle %f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f \n%c\n');
points = fscanf(fin,'triangle %f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f \n%c\n');
triangles = zeros(length(points)/16,16);
for ii = 1:length(points)/16
    triangles(ii,:) = points(16*(ii-1)+1:16*(ii-1)+16);  %x1 y1 z1 x2 y2 z2 x3 y3 z3 nx ny nz 
end

angle_x = 25;
Rx = [1 0 0;0 cosd(angle_x) -sind(angle_x);0 sind(angle_x) cosd(angle_x)];
   % iptsetpref('ImshowBorder','tight');
for yy = height-1:-1:0
    display(['trace ' num2str(yy/height*100) '%']);
    for xx = 0:width-1
        dr = zeros(1,3);
        dr(1) = (2 * ((xx+0.5)/width) -1 )*angle*aspectratio;
        dr(2) = (1 - 2*((yy+0.5)/height))*angle;
        dr(3) = 1;
        dr = dr/norm(dr);
   %     dr = (Rx*dr')';
   %     image(yy+1,xx+1,:) = raytrace([-0.1 3 0],dr,triangles,0);
        image(yy+1,xx+1,:) = raytrace([0 0 0],dr,triangles,0);
    end

    imshow(uint8(image*255));
    drawnow;
end
imshow(uint8(image*255));

