clc;
clear;
close all;
emit_photon = 100000;
estimate = 100;

global store_photonmap;
global photon_count;
angle_y = 30;
fin = fopen('prism1.asc','r');
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


      

photon_count = 0;
while (photon_count<emit_photon)
        photon_count
        
        lightn = randi([1 size(lights,1)]);
        light = lights(lightn,:);
        
        pointo = [(rand()-0.5)*light(1),(rand()-0.5)*light(2),0,1]';
        z1 = light(6:8)/norm(light(6:8));
        up_p = light(9:11) - (light(9:11)*z1')*z1;
        y1 = up_p/norm(up_p);
        x1  = cross(y1,z1);
        x_wi = [[[x1',y1',z1'],[0 0 0]'];[0 0 0 1]];
        trans = [light(1),0 0 light(3); 0, light(2),0 light(4);0 0 1 light(5);0 0 0 1]; 
        point1_t = trans*x_wi*pointo;
        origin = point1_t(1:3)/point1_t(4);
        %origin = [origin_x,origin_y,origin_z];
        nlight = light(6:8);
        
        r1 = 2*pi*rand(); r2 = rand(); r2s = sqrt(r2);
        u = cross([1 0 0],nlight);
        v = cross(nlight,u);
        direction = u*cos(r1)*r2s+v*sin(r1)*r2s+nlight*sqrt(1-r2);
        
        output_color = photonmap_caustic(origin',direction,triangles,sphere,[1 1 1],0);

        %figure(1);
        %hold on;
        %plot3(origin(1)+direction(1),origin(2)+direction(2),origin(3)+direction(3),'.');
        
%        
 

end


save('pm_caustic.mat','store_photonmap');