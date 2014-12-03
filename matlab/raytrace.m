function output_color = raytrace( origin,  direction, triangle_list, depth )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    t_min = 10000;
    triangle_index = 0;
    output_color = [0 0 0];
    for ii = 1:size(triangle_list,1)   % find the nearest intersection point
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
            for ii = 1:size(triangle_list,1)
                if (triangle_list(ii,13)~=0)   % it is a light source
                    triangle2 = triangle_list(ii,:);
                    l_position = (triangle2(1:3)+triangle2(4:6)+triangle2(7:9))/3;
                    inter_position = origin+direction*t_min;
                    dr = inter_position-l_position;
                    t_light = norm(dr);
                    dr = dr/norm(dr); 
                    T = 1;
                    
                    for jj = 1:size(triangle_list,1)
                        if (jj~=ii&& jj~=triangle_index)
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
                    output_color = output_color + max(0,-nhit*dr')*triangle(10:12).*triangle2(13:15)*T;
                end
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
                    
                    output_color = reflection*raytrace(inter_position,reflect_dir,triangle_list,depth+1)+...
                                   transmission*raytrace(inter_position,refract_dir,triangle_list,depth+1);
                    
                    
                else    % from inside to outside
                    theta1 = acosd(angle);
                    [reflection,theta2] = fresnel(n2,n1,theta1);
                    transmission = 1 - reflection;
                    reflect_dir = direction - nhit*2*(direction*nhit');
                    reflect_dir = reflect_dir/norm(reflect_dir);
                    if (transmission ==0)   % total internal reflection
                         output_color = reflection*raytrace(inter_position,reflect_dir,triangle_list,depth+1);
                    else
                        refract_dir = direction+nhit*(direction*nhit')*n1/n2-nhit*(direction*nhit');
                        refract_dir = refract_dir/norm(refract_dir);

                    
                        output_color = reflection*raytrace(inter_position,reflect_dir,triangle_list,depth+1)+...
                        transmission*raytrace(inter_position,refract_dir,triangle_list,depth+1);   
                    end
                end
            end
                
            
        end
    end
end

