function output_color = raytrace_pm( origin,  direction, triangle_list,sphere, depth ,NS, store_photonmap )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    t_min = 10000;
    triangle_index = 0;
    output_color = [0 0 0];
    
    [inter,t0,t1] = intersect_sphere(origin,direction,sphere);
    if (inter) 
        start = 1;
    else
        start = 9;
    end
    
    
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
            output_color = [0 0 0];
        elseif (triangle(16) == 68) % it is a diffusive surface
                  p1 = triangle(1:3);p2 = triangle(4:6);p3=triangle(7:9);
                  nhit = cross(p3-p2,p2-p1);
                  nhit = nhit/norm(nhit);        
                    inter_position = origin+direction*t_min;
                    
                    IDX = knnsearch(NS,inter_position,'K',500);
                    dp = sqrt(store_photonmap(IDX,1).^2+store_photonmap(IDX,2).^2+store_photonmap(IDX,3).^2);
                    
                    IDX(abs((store_photonmap(IDX,1:3)-[ones(length(IDX),1)*inter_position(1) ones(length(IDX),1)*inter_position(2) ones(length(IDX),1)*inter_position(3)])*nhit')./dp>1e-3) = [];
                    dp = sqrt(store_photonmap(IDX,1).^2+store_photonmap(IDX,2).^2+store_photonmap(IDX,3).^2);
                    k = 1.5;
                    wpc = 1-dp/k/dp(end);
                    %length(IDX)
                    color_temp = store_photonmap(IDX,4:6)'*wpc/length(IDX)/norm(store_photonmap(IDX(end),1:3)-inter_position);
                    output_color = output_color + color_temp'.*triangle(10:12);

            
        elseif (triangle(16) == 84)     % it is a transparent surface
            if (depth >= 5) 
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
                    
                    output_color = reflection*raytrace_pm(inter_position,reflect_dir,triangle_list,sphere,depth+1 ,NS, store_photonmap )+...
                                   transmission*raytrace_pm(inter_position,refract_dir,triangle_list,sphere,depth+1 ,NS, store_photonmap );
                    
                    
                else    % from inside to outside
                    theta1 = acosd(angle);
                    [reflection,theta2] = fresnel(n2,n1,theta1);
                    transmission = 1 - reflection;
                    reflect_dir = direction - nhit*2*(direction*nhit');
                    reflect_dir = reflect_dir/norm(reflect_dir);
                    if (transmission ==0)   % total internal reflection
                         output_color = reflection*raytrace_pm(inter_position,reflect_dir,triangle_list,sphere,depth+1 ,NS, store_photonmap );
                    else
                        refract_dir = direction+nhit*(direction*nhit')*n1/n2-nhit*(direction*nhit');
                        refract_dir = refract_dir/norm(refract_dir);

                    
                        output_color = reflection*raytrace_pm(inter_position,reflect_dir,triangle_list,sphere,depth+1 ,NS, store_photonmap )+...
                        transmission*raytrace_pm(inter_position,refract_dir,triangle_list,sphere,depth+1 ,NS, store_photonmap );   
                    end
                end
            end
                
            
        end
    end
end

