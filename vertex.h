#ifndef VERTEX_CLASS
#define VERTEX_CLASS

class Vertex
{
public:
	float x;
	float y;
	float z;

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

	Vertex sub(Vertex other)
	{
		Vertex v;
		v.x = x-other.x;
		v.y = y-other.y;
		v.z = z-other.z;
		return v;
	}

	Vertex add(Vertex other)
	{
		Vertex v;
		v.x = x+other.x;
		v.y = y+other.y;
		v.z = z+other.z;
		return v;
	}

	Vertex mul(Vertex other)
	{
		Vertex v;
		v.x = x*other.x;
		v.y = y*other.y;
		v.z = z*other.z;
		return v;
	}

	Vertex cross(Vertex other)
	{
		Vertex v;
		v.x = y * other.z - z*other.y;
		v.y = z * other.x - x*other.z;
		v.z = x * other.y - y*other.x;
		return v;
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
		Vertex v;
		v.x = x * s;
		v.y = y * s;
		v.z = z * s;
		return v;
	}

	void operator=(Vertex &other_vertex)
	{
		x = other_vertex.x;
		y = other_vertex.y;
		z = other_vertex.z;
	}
};

#endif //VERTEX_CLASS