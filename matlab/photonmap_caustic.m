function output_color = photonmap_caustic( origin,  direction, triangle_list, sphere,light_power,depth )

global store_photonmap;
global photon_count;

%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    t_min = 10000;
    triangle_index = 0;
    if (depth>10)
        output_color = [0 0 0];
        return;
    else
    output_color = [0 0 0];
    [inter,t0,t1] = intersect_sphere(origin,direction,sphere);
    if (depth==1 && inter == false) 
        return;
    end
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
        if (depth ==1 && triangle(16)~=84)
            return;
        end
        if (triangle(13)~=0) % it is a light source
            output_color = [0 0 0];

        elseif (triangle(16) == 68) % it is a diffusive surface

            p1 = triangle(1:3);p2 = triangle(4:6);p3=triangle(7:9);
            nhit = cross(p3-p2,p2-p1);
            nhit = nhit/norm(nhit);                 
            inter_position = origin+direction*t_min;
            
            if (depth>1)
                
                photon_count = photon_count +1;
                store_photonmap(photon_count,:) = [inter_position,light_power];
            end
            
            
           if (nhit*direction'>0)
               nhit = -nhit;
           end
           
           if (norm(light_power.*triangle(10:12))~=0)
           p = max(triangle(10:12).*light_power)/max(light_power);
           if( rand()>=p)   %absorption
               output_color = [0 0 0];
               return;
           else    % diffusive reflection
               r1 = 2*pi*rand(); r2 = rand(); r2s = sqrt(r2);
               if (abs(nhit(1))>0.1)
                   u = cross([0 1 0],nhit);
               else
                   u = cross([1 0 0],nhit);
               end
               v = cross(nhit,u);
               dir = u*cos(r1)*r2s+v*sin(r1)*r2s+nhit*sqrt(1-r2);
               
                   
                   
                output_color = photonmap_caustic(inter_position+dir*1e-2,dir,triangle_list,sphere,light_power.*triangle(10:12)/p,depth+1);
           end
           end
           
            

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
                    if (rand<transmission)
                        output_color = photonmap_caustic(inter_position,refract_dir,triangle_list,sphere,light_power*transmission,depth+1);
                    else
                        output_color = photonmap_caustic(inter_position,reflect_dir,triangle_list,sphere,light_power*reflection,depth+1);
                    end
                    
                else    % from inside to outside
                    theta1 = acosd(angle);
                    [reflection,theta2] = fresnel(n2,n1,theta1);
                    transmission = 1 - reflection;
                    reflect_dir = direction - nhit*2*(direction*nhit');
                    reflect_dir = reflect_dir/norm(reflect_dir);
                    if (transmission ==0)   % total internal reflection
                         output_color = photonmap_caustic(inter_position,reflect_dir,triangle_list,sphere,light_power,depth+1);
                    else
                        refract_dir = direction+nhit*(direction*nhit')*n1/n2-nhit*(direction*nhit');
                        refract_dir = refract_dir/norm(refract_dir);
                        if (rand<transmission)
                            output_color = photonmap_caustic(inter_position,refract_dir,triangle_list,sphere,light_power*transmission,depth+1);
                        else
                            output_color = photonmap_caustic(inter_position,reflect_dir,triangle_list,sphere,light_power*reflection,depth+1);
                        end
                    end
                end
            end
                
            
        end
     end
    
    end
end

