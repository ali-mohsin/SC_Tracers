#ifndef COLOR_CLASS
#define COLOR_CLASS

class Color 
{
public:
	float r;
	float g;
	float b;
	float a;
public:
	Color()
	{
		r = g = b = 0;
	}
	Color(float red, float green, float blue, float alpha)
	{
		r = red;
		g = green;
		b = blue;
	}
	void operator=(Color &other_color)
	{
		r = other_color.r;
		g = other_color.g;
		b = other_color.b;
	}

	Color operator+(Color &other_color)
	{
		Color c;
		c.r = r + other_color.r;
		c.g = g + other_color.g;
		c.b = b + other_color.b;

		return c;
	}
	Color mul(Color &other)
	{
		Color c;
		c.r = r*other.r;
		c.g = g*other.g;
		c.b = b*other.b;
		return c;
	}
	Color scaleColor(float s)
	{
		Color c;
		c.r = r * s;
		c.g = g * s;
		c.b = b * s;
		return c;
	}
	
	short ctoi(float color)		/* convert float color to GzIntensity short */
	{
		return(short)((int)(color * ((1 << 12) - 1)));
	}
};

#endif //COLOR_CLASS