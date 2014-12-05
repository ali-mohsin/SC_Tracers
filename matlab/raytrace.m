function output_color = raytrace( origin,  direction, triangle_list,sphere,lights, depth )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    
x_l = linspace(-0.5,0.5,5);
y_l = linspace(-0.5,0.5,5);
%x_l = 0;
%y_l = 0;
count = 0;
light_list = zeros(length(x_l)*length(y_l)*size(lights,1),3);
for l_ii = 1:size(lights,1);
    light = lights(l_ii,:);

    trans = [light(1),0 0 light(3); 0, light(2),0 light(4);0 0 1 light(5);0 0 0 1]; 

    z1 = light(6:8)/norm(light(6:8));
    up_p = light(9:11) - (light(9:11)*z1')*z1;
    y1 = up_p/norm(up_p);
    x1  = cross(y1,z1);
    x_wi = [[[x1',y1',z1'],[0 0 0]'];[0 0 0 1]];
    

    for ii = 1:length(x_l)
        for jj = 1:length(y_l)
            pos = [x_l(ii), y_l(jj) ,0, 1]';
            pos1 = trans*x_wi*pos;
            xxx = pos1(1)/pos1(4);
            yyy = pos1(2)/pos1(4);
            zzz = pos1(3)/pos1(4);
            count = count+1;
            light_list(count,:) = [xxx,yyy,zzz];
       
        end
    end
end
[inter,t0,t1] = intersect_sphere(origin,direction,sphere);
if (inter) 
    start = 1;
else
    start = 9;
end
t_min = 10000;
    triangle_index = 0;
    output_color = [0 0 0];
    for ii = start:size(triangle_list,1)   % find the nearest intersection point
        triangle = triangle_list(ii,:);
        [inter,t] = intersect_tri(origin,direction,triangle);

        if ((inter==true) && (t<t_min) && (t>1e-3)) 
            
            t_min = t;
            triangle_index = ii;
        end
    end
    
    
    if (triangle_index~=0)   % if some triangle is found intersected by the ray
        triangle = triangle_list(triangle_index,:);
        
        if (triangle(13)~=0) % it is a light source
            output_color = triangle(10:12);
        elseif (triangle(16) == 68) % it is a diffusive surface
            
            for ii = 1:size(light_list,1)
               
                    l_position = light_list(ii,:);
                    
                    inter_position = origin+direction*t_min;
                    dr = inter_position-l_position;
                    t_light = norm(dr);
                    dr = dr/norm(dr); 
                    T = 1;
                    [inter,t0,t1] = intersect_sphere(inter_position,-dr,sphere);
                    if (inter)
                        start = 1;
                    else
                        start = 9;
                    end
                    
                    for jj = start:size(triangle_list,1)
                        if (jj~=triangle_index && triangle_list(jj,13)==0)
                            triangle_test = triangle_list(jj,:);
                            [inter,t] = intersect_tri(inter_position,-dr,triangle_test);
                            if (inter==true && t<t_light && t>1e-3)
                                T = 0;
                                break;
                            end
                        end
                    end
                    p1 = triangle(1:3);p2 = triangle(4:6);p3=triangle(7:9);
                    nhit = cross(p3-p2,p2-p1);
                    nhit = nhit/norm(nhit);                  
                    output_color = output_color + max(0,-nhit*dr')*triangle(10:12).*[1 1 1]*T/size(light_list,1);
                
            end
            
        elseif (triangle(16) == 84)     % it is a transparent surface
            if (depth > 10) 
                return;
            else
                p1 = triangle(1:3);p2 = triangle(4:6);p3=triangle(7:9);
                nhit = cross(p3-p2,p2-p1);
                nhit = nhit/norm(nhit); 
                n1 = 1.0;
                n2 = 2.4;
                % n2= 2.43;
               %  n2 = 2.46;
                angle = direction*nhit';
                inter_position = origin+direction*t_min;
                
                if (angle<0)   % from outside to inside
                    theta1 = acosd(-angle);
                    [reflection,theta2] = fresnel(n1,n2,theta1);
                    transmission = 1 - reflection;
                  %  direction
                    reflect_dir = direction - nhit*2*(direction*nhit');
                    reflect_dir = reflect_dir/norm(reflect_dir);
                    refract_dir = direction+nhit*(direction*nhit')*n2/n1-nhit*(direction*nhit');
                    refract_dir = refract_dir/norm(refract_dir);
                   % refract_dir
                  %  depth
                    
                    output_color = reflection*raytrace(inter_position,reflect_dir,triangle_list,sphere,lights,depth+1)+...
                                   transmission*raytrace(inter_position,refract_dir,triangle_list,sphere,lights,depth+1);
                    
                    
                else    % from inside to outside
                    theta1 = acosd(angle);
                    [reflection,theta2] = fresnel(n2,n1,theta1);
                    transmission = 1 - reflection;
                    reflect_dir = direction - nhit*2*(direction*nhit');
                    reflect_dir = reflect_dir/norm(reflect_dir);
                    if (transmission ==0)   % total internal reflection
                         output_color = reflection*raytrace(inter_position,reflect_dir,triangle_list,sphere,lights,depth+1);
                    else
                        refract_dir = direction+nhit*(direction*nhit')*n1/n2-nhit*(direction*nhit');
                        refract_dir = refract_dir/norm(refract_dir);

                    
                        output_color = reflection*raytrace(inter_position,reflect_dir,triangle_list,sphere,lights,depth+1)+...
                        transmission*raytrace(inter_position,refract_dir,triangle_list,sphere,lights,depth+1);   
                    end
                end
            end
                
            
        end
    end
end

