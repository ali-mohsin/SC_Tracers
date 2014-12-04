#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <assert.h>

using namespace std;

#include "fresnel.h"
#include "parser.h"
#include "diamond.h"
#include "raytracer.h"

#define WIDTH 256
#define HEIGHT 256
#define FOV  20.0
#define ASPECTRATIO WIDTH/HEIGHT
#define ANGLE tand(0.5*FOV)

#define FAIL 0
#define SUCCESS 1


#define	AAKERNEL_SIZE	6

#define Normalize255(x) min( max(short(0),x), short(4095) )>>4 

short	ctoi(float color)	
{
	return (short)((int)(color * ((1 << 12) - 1)));
}

typedef	struct {
	float    red;	
	float    green;
	float    blue;
} pixel;

class ImageBuffer{
public:
	int width;
	int height;
	pixel *buf;
public:
	ImageBuffer( int width, int height)
	{
		this->width = width;
		this->height = height;
		buf = new pixel[ width * height ];
	}

	int initImageBuf()
	{
		for( short y = 0 ; y < height ; y++ )
			for( short x = 0 ; x < width ; x++ )
			{
				int pixelNum = x + y*width;
				buf[pixelNum].red = 0;
				buf[pixelNum].green = 0;
				buf[pixelNum].blue = 0;
			}
			return SUCCESS;
	}


	int setImageBuf(int i, int j, float r, float g, float b)
	{
		/* write pixel values into the image buffer */

		i = max(0, min(i, width));
		j = max(0, min(j, height));
		int pixelNum = i + j*width; 

		buf[pixelNum].red = (r);
		buf[pixelNum].green =(g);
		buf[pixelNum].blue = (b);

		if(buf[pixelNum].red > 1)
			buf[pixelNum].red = 1;
		if(buf[pixelNum].green > 1)
			buf[pixelNum].green = 1;
		if(buf[pixelNum].blue > 1)
			buf[pixelNum].blue = 1;
		if(buf[pixelNum].red <0)
			buf[pixelNum].red = 0;
		if(buf[pixelNum].green<0)
			buf[pixelNum].green = 0;
		if(buf[pixelNum].blue <0)
			buf[pixelNum].blue = 0;



		return SUCCESS;
	}

	int addImageBuf(int i, int j, float r, float g, float b)
	{
		/* write pixel values into the image buffer */

		i = max(0, min(i, width));
		j = max(0, min(j, height));
		int pixelNum = i + j*width; 

		buf[pixelNum].red += (r);
		buf[pixelNum].green +=(g);
		buf[pixelNum].blue += (b);

		if(buf[pixelNum].red > 1)
			buf[pixelNum].red = 1;
		if(buf[pixelNum].green > 1)
			buf[pixelNum].green = 1;
		if(buf[pixelNum].blue > 1)
			buf[pixelNum].blue = 1;
		if(buf[pixelNum].red <0)
			buf[pixelNum].red = 0;
		if(buf[pixelNum].green<0)
			buf[pixelNum].green = 0;
		if(buf[pixelNum].blue <0)
			buf[pixelNum].blue = 0;



		return SUCCESS;
	}


	int flushDisplay2File(FILE* outfile)
	{
		int i, j;
		(void) fprintf(outfile, "P6\n%d %d\n255\n", WIDTH, HEIGHT);
		for (j = 0; j < WIDTH; ++j)
		{
			for (i = 0; i < HEIGHT; ++i)
			{
				int pixelNum = i+j*width;
				static unsigned char color[3];
				color[0] = buf[pixelNum].red *255;  /* red */
				color[1] = buf[pixelNum].green * 255;  /* green */
				color[2] = buf[pixelNum].blue *255;  /* blue */
				(void) fwrite(color, 1, 3, outfile);
			}
		}
		(void) fclose(outfile);
		return 1;
	}

	int freeBuf()
	{
		delete buf;
		return SUCCESS;
	}
};

void matMult(float m1[][4],float m2[][4],float ans[][4])
	{
		for(int i=0;i<4;i++)
		{
			for(int j=0;j<4;j++)
			{
				ans[i][j]=0;
			}
		}

		for(int i=0;i<4;i++)
		{
			for(int j=0;j<4;j++)
			{
				for(int k=0;k<4;k++)
				{
					ans[i][j] += m1[i][k] * m2[k][j]; 

				}
			}
		}
	}


void pointMult(float m1[][4],float v[4],float ans[4])
	{
		float aans[4];
		for(int i=0;i<4;i++)
		{
				aans[i]=0;
		}

		for(int i=0;i<4;i++)
		{
			for(int j=0;j<4;j++)
			{
				aans[i] += m1[i][j] * v[j]; 
			}
		}

		for(int i=0;i<4;i++)
			ans[i]=aans[i];
	}


int main ( int argc)
{
	Triangle t1;
	Triangle t2;
	int	status = 0; 

	int lights = 4; //no.of lights
	ImageBuffer image(WIDTH,HEIGHT);
	status |= image.initImageBuf();
	Parser* p = new Parser("E:\\USC\\CSCI 580 (Fall 2014)\\Final\\diamond\\diamond_trans_330.asc");
	Diamond d = Diamond(p->parse());
	int* count= new int(0);
	float angle_x = 25;
	float Rx[3][3] =
    {
		{1,0,0},
		{0,cosd(angle_x),-sind(angle_x)},
		{0,sind(angle_x),cosd(angle_x)}
    };


	float AAFilter[AAKERNEL_SIZE][3] 	=		/* X, Y, coef */
	{
		-0.52, 0.38, 0.128,
		0.41, 0.56, 0.119,
		0.27, 0.08, 0.294,
		-0.17, -0.29, 0.249,
		0.58, -0.55, 0.104,
		-0.31, -0.71, 0.106
	};


	Vertex dir = Vertex(0,0,10.75);//-[-0.1,3,6]);
	Vertex temp= Vertex(-0.1,3,6);
	dir=dir.sub(temp);
	dir.norm();// = direction/norm(direction);
	const int num_lights=2;

	float light_vec[][11] = 
	{
		{1,1,0,2.745,10.75, 0, -1, 0, 0, 0, -1},
		{1,1,0,0,6,dir.x,dir.y,dir.z,0,1,0}
	};


	float x_l[5] = {-0.5, -0.25, 0, 0.25, 0.5};
	float y_l[5] = {-0.5, -0.25, 0, 0.25, 0.5};
	float light_list[5*5*num_lights][3] = {0}; // zeros(num_lights,3);

	

	int count_light =0;
	for(int i=0;i<num_lights;i++)
	{
		float point1[4] = {-0.5,-0.5,0,1};//'';
	    float point2[4] ={0.5,-0.5,0,1};//]';
        float point3[4] = {0.5,0.5,0,1};//]';
        float point4[4] = {-0.5,0.5,0,1};//]';


		float cur_light[11] ;
		for(int k=0;k<11;k++)
			cur_light[k]=light_vec[i][k];
    
		Vertex z1 = Vertex(cur_light[5],cur_light[6],cur_light[7]);///norm(light(6:8));
		z1.norm();

		Vertex up_p = Vertex(cur_light[8],cur_light[9],cur_light[10]);///norm(light(6:8));
		up_p=up_p.mul(z1);
		up_p=up_p.mul(z1);
		Vertex y1=Vertex(cur_light[8],cur_light[9],cur_light[10]);
		y1=y1.sub(up_p);
		y1.norm();
		Vertex x1=y1.cross(z1);

		float x_wi[4][4]=
		{
			{x1.x,y1.x,z1.x,0},
			{x1.y,y1.y,z1.y,0},
			{x1.z,y1.z,z1.z,0},
			{0,0,0,1}
		};

//    x_wi = [[[x1',y1',z1'],[0 0 0]'];[0 0 0 1]];
		float trans[4][4] = 
		{
			{light_vec[i][0],0, 0 ,light_vec[i][2]},
			{ 0, light_vec[i][1],0 ,light_vec[i][3]},
			{0 ,0, 1, light_vec[i][4]},
			{0,0,0,1}
		}; 

		float M[4][4];

		matMult(trans,x_wi,M);
		pointMult(M,point1,point1);
		pointMult(M,point2,point2);
		pointMult(M,point3,point3);
		pointMult(M,point4,point4);
    
		for(int x=0;x<3;x++)
		{
			point1[x]/=point1[3];
			point2[x]/=point2[3];
			point3[x]/=point3[3];
			point4[x]/=point4[3];
		}
    
		Vertex vv1 = Vertex(point1[0],point1[1],point1[2]);
		Vertex vv2 = Vertex(point2[0],point2[1],point2[2]);
		Vertex vv3 = Vertex(point3[0],point3[1],point3[2]);
		Vertex vv4 = Vertex(point4[0],point4[1],point4[2]);

		 t1 = Triangle (vv1,vv2,vv4,"L");
		 t2 = Triangle (vv4,vv2,vv3,"L");
		t1.surfaceColor = Color(1,1,1,1);
		t1.emissionColor = Color(1,1,1,1);
		t2.surfaceColor = Color(1,1,1,1);
		t2.emissionColor = Color(1,1,1,1);
		d.tri.push_back(t1);
		d.tri.push_back(t2);

		for (int ii=0;ii<5;ii++){
			for(int jj=0;jj<5;jj++){
				float pos[4]={x_l[ii],y_l[jj],0,1};	
				float M[4][4];
				matMult(trans,x_wi,M);
				pointMult(M,pos,pos);
				pos[0]/=pos[3];
				pos[1]/=pos[3];
				pos[2]/=pos[3];

				light_list[count_light][0]= pos[0];
				light_list[count_light][1]= pos[1];
				light_list[count_light][2]= pos[2];
				count_light++;
			}
		}


	}

/*
	float x_l[5] =  {-0.5000,   -0.2500  ,       0   , 0.2500  ,  0.5000};
	float y_l[5] =  {-0.5000,   -0.2500  ,       0   , 0.2500  ,  0.5000};
	int count = 0;t
	float light_list[5*5*num_lights][3] = {0};// zeros(num_lights,3);
	for (int li=0; li < num_lights;li++)
	{
		
	}
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
*/

	for ( int y = 0 ; y < HEIGHT ; y++ )
	{
		//display(['trace ' num2str(yy/height*100) '%']);
		cout<<"tracing "<< (double)y/(double)HEIGHT*100 <<"%"<<endl;
		for ( int x = 0 ; x < WIDTH ; x++ )
		{
			for (int aa = 0; aa < AAKERNEL_SIZE; aa++)
			{
				Vertex dr(0,0,0);
				dr.x = (2 * ((x+0.5 + AAFilter[aa][0])/WIDTH) -1  )*ANGLE*ASPECTRATIO;
				float temp = (1 - 2*((y+0.5 + AAFilter[aa][1])/HEIGHT)) ;
				float angle=tand(0.5*FOV);;
				dr.y=temp*angle;
				dr.z = 1;
				float n=normalize(dr);
				dr.x = dr.x/n;	dr.y = dr.y/n;	dr.z = dr.z/n;
				Vertex v(0,0,0);
				Color output_color_r,output_color_g,output_color_b;
				dr.matMul(Rx);
				Vertex dr_origin;
				dr_origin.x=dr.x*9-0.1;
				dr_origin.y=dr.y*9+3;
				dr_origin.z=dr.z*9;
				output_color_r = raytrace(dr_origin,dr,d.tri,0,count,count_light,2.40,light_list);
				output_color_g = raytrace(dr_origin,dr,d.tri,0,count,count_light,2.43,light_list);
				output_color_b = raytrace(dr_origin,dr,d.tri,0,count,count_light,2.46,light_list);

				image.addImageBuf(x,y, (output_color_r.r)*AAFilter[aa][2], (output_color_g.g)*AAFilter[aa][2], (output_color_b.b)*AAFilter[aa][2]);
				//image.addImageBuf(x,y, (output_color_r.r)*AAFilter[aa][2], (output_color_r.g)*AAFilter[aa][2], (output_color_r.b)*AAFilter[aa][2]);
		
			}
		}
	}
	*count=*count+1-1;

	FILE *outfile;
	if( (outfile  = fopen( "330.ppm" , "wb" )) == NULL )
	{
		cout<<"Failed to open output file"<<endl;
		return 0;
	}

	image.flushDisplay2File(outfile);
	image.freeBuf();
}