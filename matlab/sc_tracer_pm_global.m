clc;
clear;
close all;
angle_y = 30;


width = 256;
height = 256;

image = zeros(height,width,3);

FOV = 40;
aspectratio = width/height;

AAFilter=...
	[...
		-0.52, 0.38, 0.128;...
		0.41, 0.56, 0.119;...
		0.27, 0.08, 0.294;...
		-0.17, -0.29, 0.249;...
		0.58, -0.55, 0.104;...
		-0.31, -0.71, 0.106...
	];

angle = tand(0.5*FOV);
load('pm_global.mat');
NS = KDTreeSearcher(store_photonmap(:,1:3),'Distance','euclidean');
fin = fopen('prism.asc','r');
points = fscanf(fin,'triangle %f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f \n%c\n');
count = 0;
triangles = zeros(length(points)/16,16);
for ii = 1:length(points)/16
    count = count+1;
    triangles(count,:) = points(16*(ii-1)+1:16*(ii-1)+16);  %x1 y1 z1 x2 y2 z2 x3 y3 z3 nx ny nz 
end

fclose(fin);

sphere = [0,-0.75,10.75,2.1];
fin = fopen('room.asc','r');
points = fscanf(fin,'triangle %f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f \n%c\n');

for ii = 1:length(points)/16
    count = count+1;
    triangles(count,:) = points(16*(ii-1)+1:16*(ii-1)+16);  %x1 y1 z1 x2 y2 z2 x3 y3 z3 nx ny nz 
end
fclose(fin);



% width, length, position, direction, world up
direction = ([0,0,10.75]-[-0.1,3,6]);
direction = direction/norm(direction);
lights = [1.2,1.2,1.6,2.745,10.75, 0, -1, 0, 0, 0, -1];
%%
for ii = 1:size(lights,1)
    point1 = [-0.5,-0.5,0,1]';
    point2 = [0.5,-0.5,0,1]';
    point3 = [0.5,0.5,0,1]';
    point4 = [-0.5,0.5,0,1]';
    light = lights(ii,:);
    z1 = light(6:8)/norm(light(6:8));
    up_p = light(9:11) - (light(9:11)*z1')*z1;
    y1 = up_p/norm(up_p);
    x1  = cross(y1,z1);
    x_wi = [[[x1',y1',z1'],[0 0 0]'];[0 0 0 1]];
    trans = [light(1),0 0 light(3); 0, light(2),0 light(4);0 0 1 light(5);0 0 0 1]; 
    point1_t = trans*x_wi*point1;
    point2_t = trans*x_wi*point2;
    point3_t = trans*x_wi*point3;
    point4_t = trans*x_wi*point4;
    
    point1_t(1:3) = point1_t(1:3)/point1_t(4);
    point2_t(1:3) = point2_t(1:3)/point2_t(4);
    point3_t(1:3) = point3_t(1:3)/point3_t(4);
    point4_t(1:3) = point4_t(1:3)/point4_t(4);
    
    
    
    triangles_temp1 = [point1_t(1:3)',point2_t(1:3)',point4_t(1:3)',1,1,1,1,1,1,68];
    triangles_temp2 = [point4_t(1:3)',point2_t(1:3)',point3_t(1:3)',1,1,1,1,1,1,68];
    triangles = [triangles;triangles_temp1;triangles_temp2];
end

%%
angle_x = 25;
Rx = [1 0 0;0 cosd(angle_x) -sind(angle_x);0 sind(angle_x) cosd(angle_x)];
for yy = height-1:-1:0
    display(['trace ' num2str(yy/height*100) '%']);
    for xx = 0:width-1
        for aa = 1:6
        dr = zeros(1,3);
        dr(1) = (2 * ((xx+0.5)/width) -1 )*angle*aspectratio;
        dr(2) = (1 - 2*((yy+0.5)/height))*angle;
        dr(3) = 1;
        dr = dr/norm(dr);
        %dr = (Rx*dr')';
        color_temp = raytrace_pm([0 0 0],dr,triangles,sphere,0,NS,store_photonmap)*AAFilter(aa,3);
        image(yy+1,xx+1,1) = image(yy+1,xx+1,1)+ color_temp(1);
        image(yy+1,xx+1,2) = image(yy+1,xx+1,2)+ color_temp(2);
        image(yy+1,xx+1,3) = image(yy+1,xx+1,3)+ color_temp(3);
        %imshow(uint8(image*255));
        %drawnow;
    
        end
    end
            imshow(uint8(image*255));
        drawnow;
end
%imshow(uint8(image*255));

save(['image_global.mat'],'image');
