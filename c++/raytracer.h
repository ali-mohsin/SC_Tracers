#include "triangle.h"
#include "vertex.h"
#include "math.h"

//#define Normalize255(x) min( max(short(0),x), short(4095) )>>4 
//
//short	ctoi(float color)	
//{
//	return(short)((int)(color * ((1 << 12) - 1)));
//}
float normalize(Vertex v)
{
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

Color raytrace(Vertex origin, Vertex direction, vector<Triangle> triangle_list, float depth, int* count,int lights)
{
	float t_min = 10000;
	int triangle_index = -1;
	Color output_color(0,0,0,1);
	Triangle tri1, tri2;
	for( int i = 0; i < triangle_list.size(); i++ )
	{
		tri1 = triangle_list[i];
		tri1.intersection(origin, direction);
		bool inter = tri1.inter;
		float t = tri1.t;
		if( (inter == true) && (t < t_min) && (t>0.001) )
		{
			t_min = t;
			triangle_index = i;
		}
	}

	if( triangle_index != -1)			// if some triangle is found intersected by the ray
	{
		*count=*count+1;
		tri1 = triangle_list[triangle_index];
		if ( tri1.token == "L" )
		{
			output_color = tri1.surfaceColor;
		//	output_color.r *= 255; 
		//	output_color.g *= 255;
		//	output_color.b *= 255;
			return output_color;
		}
		else  if(tri1.token == "D")
			{ 
				// it is a diffusive surface
			for( int i = 0; i < triangle_list.size(); i++ )
			{
				if( triangle_list[i].token == "L" )		// it is a light surface
				{
					tri2 = triangle_list[i];
					Vertex l_position(0,0,0), inter_position(0,0,0);
					l_position += tri2.v1;
					l_position += tri2.v2;
					l_position += tri2.v3;
					l_position.x /= 3;	l_position.y /= 3;	l_position.z /= 3;
					inter_position.x = origin.x + direction.x*t_min;	inter_position.y = origin.y + direction.y*t_min;	inter_position.z = origin.z + direction.z*t_min;
					Vertex dr = inter_position.sub(l_position);
					float t_light = normalize(dr);
					dr.x = dr.x/t_light;  dr.y = dr.y/t_light;  dr.z = dr.z/t_light;			// dr = dr/norm(dr)
					float T = (float)1/lights; //normalizing so that we divide by number of light sources.
					for ( int j = 0; j < triangle_list.size(); j++ )
					{
                        if ( j != i && j != triangle_index )
						{
                            Triangle triangle_test = triangle_list[j];
							Vertex neg_dr(-dr.x,-dr.y,-dr.z);
							triangle_test.intersection(inter_position, neg_dr);
							bool inter = triangle_test.inter;
							float t = triangle_test.t;
                            if (inter == true && t < t_light && t > 0.001)
							{
                                T = 0;
                                break;
							}
						}
					}
                    Vertex p1 = tri1.v1 , p2 = tri1.v2 , p3 = tri1.v3;
					p3 -= p2;
					p2 -= p1;
                    Vertex nhit = p3.cross(p2);		// cross(p3-p2,p2-p1)
					float nhit_norm = normalize(nhit);
                    nhit.x = nhit.x/nhit_norm;     nhit.y = nhit.y/nhit_norm;   nhit.z = nhit.z/nhit_norm;             // nhit = nhit/norm(nhit) 
					Vertex neg_nhit(-nhit.x,-nhit.y,-nhit.z);
                    output_color = output_color + tri1.surfaceColor.scaleColor(max(float(0),neg_nhit.dot(dr))).mul(tri2.emissionColor.scaleColor(T));		
				}
			}
			return output_color;
		}
        else if (tri1.token == "T")			// It is a transparent surface
		{
            if (depth > 5)
                return Color(0,0,0,1);
            else
			{
                Vertex p1 = tri1.v1 , p2 = tri1.v2 , p3 = tri1.v3;
                p3 -= p2;
                p2 -= p1;
                Vertex nhit = p3.cross(p2);
                float nhit_norm = normalize(nhit);
                nhit.x = nhit.x/nhit_norm;    nhit.y = nhit.y/nhit_norm;  nhit.z = nhit.z/nhit_norm;             // nhit = nhit/norm(nhit)   

                float n1 = 1.0;
                float n2 = 2.4;

                float angle = nhit.dot(direction);
                Vertex inter_position(0,0,0);
                inter_position.x = origin.x + direction.x*t_min;  inter_position.y = origin.y + direction.y*t_min;  inter_position.z = origin.z + direction.z*t_min;

                if ( angle < 0 )		// from outside to inside
				{
                    float theta1 = acosd(-angle);
                    Fresnel output = fresnel(n1,n2,theta1);
                    float reflection = output.reflection;
                    float theta2 = output.theta2;
                    float transmission = 1 - reflection;
                    Vertex reflect_dir = direction.sub((nhit.scaleVertex(2)).scaleVertex(direction.dot(nhit))); 
                    float reflect_norm = normalize(reflect_dir);
                    reflect_dir.x = reflect_dir.x/reflect_norm;
                    reflect_dir.y = reflect_dir.y/reflect_norm;
                    reflect_dir.z = reflect_dir.z/reflect_norm;

                    Vertex refract_dir = (direction.add((nhit.scaleVertex(direction.dot(nhit))).scaleVertex(n2/n1))).sub(nhit.scaleVertex(direction.dot(nhit)));
                    float refract_norm = normalize(refract_dir);
                    refract_dir.x = refract_dir.x/refract_norm;
                    refract_dir.y = refract_dir.y/refract_norm;
                    refract_dir.z = refract_dir.z/refract_norm;
                    output_color = raytrace(inter_position,reflect_dir,triangle_list,depth+1,count,lights).scaleColor(reflection) +
                                   raytrace(inter_position,refract_dir,triangle_list,depth+1,count,lights).scaleColor(transmission);
					return output_color;
                }
                else			// from inside to outside
				{
                    float theta1 = acosd(angle);
                    Fresnel output = fresnel(n2,n1,theta1);
                    float reflection = output.reflection;
                    float theta2 = output.theta2;
                    float transmission = 1- reflection;
                    Vertex reflect_dir = direction.sub(nhit.scaleVertex(2).scaleVertex(direction.dot(nhit))); 
                    float reflect_norm = normalize(reflect_dir);
                    reflect_dir.x = reflect_dir.x/reflect_norm;
                    reflect_dir.y = reflect_dir.y/reflect_norm;
                    reflect_dir.z = reflect_dir.z/reflect_norm;
                    if (transmission == 0)			//total internal reflection
					{
                        output_color = raytrace(inter_position,reflect_dir,triangle_list,depth+1,count,lights).scaleColor(reflection);
						return output_color;
					}
                    else
                    {
						Vertex refract_dir = (direction.add((nhit.scaleVertex(direction.dot(nhit))).scaleVertex(n1/n2))).sub(nhit.scaleVertex(direction.dot(nhit)));
                        float refract_norm = normalize(refract_dir);
                        refract_dir.x = refract_dir.x/refract_norm;
                        refract_dir.y = refract_dir.y/refract_norm;
                        refract_dir.z = refract_dir.z/refract_norm;
                        output_color = raytrace(inter_position,reflect_dir,triangle_list,depth+1,count,lights).scaleColor(reflection) + raytrace(inter_position,refract_dir,triangle_list,depth+1,count,lights).scaleColor(transmission);
						return output_color;
					}
                }
			}
        }
	}
	return Color(0,0,0,1);
}
    
    
