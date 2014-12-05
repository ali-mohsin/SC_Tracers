function  [inter,t0,t1] = intersect_sphere (origin,direction,sphere)
    center = sphere(1:3);
    radius = sphere(4);
    radius2 = radius*radius;
    l = center - origin;
    tca = l*direction';
    d2 = norm(l)^2-tca^2;
    if (d2 > radius2) 
        inter = false;
        t0 = 0;
        t1 = 0;
        return;
    end
    thc = sqrt(radius2-d2);
    t0 = tca - thc;
    t1 = tca + thc;
    inter =true;
end