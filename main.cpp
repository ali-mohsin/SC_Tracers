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

#include "parser.h"
#include "diamond.h"
#include "raytracer.h"

#define WIDTH 256.0
#define HEIGHT 256.0
#define FOV  40.0
#define ASPECTRATIO WIDTH/HEIGHT
#define ANGLE tand(0.5*FOV)

#define FAIL 0
#define SUCCESS 1

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

		buf[pixelNum].red = (r);
		buf[pixelNum].green =(g);
		buf[pixelNum].blue = (b);

		return SUCCESS;
	}

	int flushDisplay2File(FILE* outfile)
	{
		int i, j;
		(void) fprintf(outfile, "P6\n%d %d\n255\n", 256, 256);
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

int main()
{
	int	status = 0; 

	int lights = 8; //no.of lights

	ImageBuffer image(WIDTH,HEIGHT);
	status |= image.initImageBuf();
	Parser* p = new Parser("prism3.asc");
	Diamond d = Diamond(p->parse());
	int* count= new int(0);
	for ( int y = 0 ; y < 256 ; y++ )
	{
		//display(['trace ' num2str(yy/height*100) '%']);
		for ( int x = 0 ; x < 256 ; x++ )
		{
			Vertex dr(0,0,0);
			dr.x = (2 * ((x+0.5)/WIDTH) -1 )*ANGLE*ASPECTRATIO;
			float temp = (1 - 2*((y+0.5)/HEIGHT));
			float angle=tand(0.5*FOV);;
			dr.y=temp*angle;
			dr.z = 1;
			float n=normalize(dr);
			dr.x = dr.x/normalize(dr);	dr.y = dr.y/normalize(dr);	dr.z = dr.z/normalize(dr);
			Vertex v(0,0,0);
			Color output_color;
			output_color = raytrace(v,dr,d.tri,0,count,lights);



			image.setImageBuf(x,y, (output_color.r), (output_color.g), (output_color.b));
		}
	}
	*count=*count+1-1;

	FILE *outfile;
	if( (outfile  = fopen( "DIAMOND.ppm" , "wb" )) == NULL )
	{
		cout<<"Failed to open output file"<<endl;
		return 0;
	}

	image.flushDisplay2File(outfile);
	image.freeBuf();
}