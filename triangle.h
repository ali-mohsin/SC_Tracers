#ifndef TRIANGLE_CLASS
#define TRIANGLE_CLASS

#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <assert.h>

#include "color.h"
#include "vertex.h"

using namespace std;

class Triangle
{
public:	
	Vertex v1;
	Vertex v2;
	Vertex v3;
	Vertex edge1;
	Vertex edge2;
	Vertex norm;
	float t; //output of intersection
	bool inter; // same
	Color surfaceColor;
	Color emissionColor;
	string token;		// T or D or L- transparent or diffusive surface or light source surface

	Triangle()
	{
		v1 = Vertex(0,0,0);
		v2 = Vertex(0,0,0);
		v3 = Vertex(0,0,0);
		token = "T";
		surfaceColor = Color(0,0,0,1);
		emissionColor =  Color(0,0,0,1);
		computeNormal();
	}

	Triangle(Vertex a,Vertex b,Vertex c, string tld_token)
	{
		v1 = a;
		v2 = b;
		v3 = c;
		token = tld_token;
		surfaceColor = Color(0,0,0,1);
		emissionColor =  Color(0,0,0,1);
		computeNormal();
	}

	void computeNormal()
	{
		edge1 = v2.sub(v1);
		edge2 = v3.sub(v1);
		norm = edge1.cross(edge2);
	}

	void print()
	{
		cout<<"Vertex 1: "<<v1.x<<" "<<v1.y<<" "<<v1.z<<endl;
		cout<<"Vertex 2: "<<v2.x<<" "<<v2.y<<" "<<v2.z<<endl;
		cout<<"Vertex 3: "<<v3.x<<" "<<v3.y<<" "<<v3.z<<endl;
		cout<<"Normal : "<<norm.x<<" "<<norm.y<<" "<<norm.z<<endl;

	}

	void intersection(Vertex origin, Vertex dir)
	{
		Vertex pvec = dir.cross(edge2);
		float det=edge1.dot(pvec);
		if(det==0)
		{
			inter=false;
			t=0;
			return;
		}

		float invDet=1/det;
		Vertex tvec=origin.sub(v1);
		float u=(tvec.dot(pvec))*(invDet);

		if(u<0 || u>1)
		{
			inter=false;
			t=0;
			return;
		}

		Vertex qvec = tvec.cross(edge1);
		float v = dir.dot(qvec)*(invDet);

		if (v<0||v+u>1)
        {
        	inter = false; 
        	t = 0; 
        	return;
        }

        inter = true;
        t = (edge2.dot(qvec))*(invDet);
	}

	void setColor(string type, Color color)
	{
		if ( type == "surface" )
			surfaceColor = color;
		else if ( type == "emission" )
			emissionColor = color;
	}
};

#endif //TRIANGLE_CLASS