#ifndef VERTEX_CLASS
#define VERTEX_CLASS

class Vertex
{
public:
	float x;
	float y;
	float z;
	float v[3];

	Vertex(float a,float b,float c)
	{
		x=a;
		y=b;
		z=c;

	}

	Vertex()
	{
		x = y = z = 0;
	}

	void operator-=(Vertex &other)
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;
	}
	void operator+=(Vertex &other)
	{
		x += other.x;
		y += other.y;
		z += other.z;
	}

	void operator=(Vertex &other_color)
	{
		x = other_color.x;
		y = other_color.y;
		z = other_color.z;
	}

	Vertex sub(Vertex other)
	{
		Vertex vv;
		vv.x = x-other.x;
		vv.y = y-other.y;
		vv.z = z-other.z;
		return vv;
	}

	Vertex add(Vertex other)
	{
		Vertex vv;
		vv.x = x+other.x;
		vv.y = y+other.y;
		vv.z = z+other.z;
		return vv;
	}

	Vertex mul(Vertex other)
	{
		Vertex vv;
		vv.x = x*other.x;
		vv.y = y*other.y;
		vv.z = z*other.z;
		return vv;
	}

	Vertex cross(Vertex other)
	{
		Vertex vv;
		vv.x = y * other.z - z*other.y;
		vv.y = z * other.x - x*other.z;
		vv.z = x * other.y - y*other.x;
		return vv;
	}

	float dot(Vertex other)
	{
		return x*other.x+y*other.y+z*other.z;
	}

	float scale(float s)
	{
		return x*s+y*s+z*s;
	}

	Vertex scaleVertex(float s)
	{
		Vertex vv;
		vv.x = x * s;
		vv.y = y * s;
		vv.z = z * s;
		return vv;
	}



	void matMul(float m[][3])
	{
		v[0]=x;
		v[1]=y;
		v[2]=z;
		float ans[3];
		for(int i=0;i<3;i++)
		{
			ans[i]=0;
		}

		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				//std::cout<<m[i][j]<<" "<<v[j]<<"\n";
				ans[i]+=m[i][j]*v[j];
			}
		}

		x=ans[0];
		y=ans[1];
		z=ans[2];
	}

	void norm()
	{
		float mag=sqrt(x*x+y*y+z*z);
		x/=mag;
		y/=mag;
		z/=mag;
	
	}
};

#endif //VERTEX_CLASS