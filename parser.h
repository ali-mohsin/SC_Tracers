#include <fstream>
#include "triangle.h"
#include "vertex.h"
#include "color.h"

class Parser
{
public:
	string file_name;
	vector<Triangle> triangles;
	ifstream fin;

	Parser(string f)
	{
		file_name = f;
	}

	vector<Triangle> parse()
	{
		fin.open(file_name);
		cout<<"----------------"<<endl;
		if(fin.fail())
		{
			cout<<"Could not open file"<<endl;
			exit(1);
		}

		string buffer;
	    vector<float> v;
		while(!fin.eof())
		{
			getline(fin,buffer);
		    istringstream buf(buffer);
		    for(string token; getline(buf, token,' '); )
	        {
	        	if (token=="triangle")
	        		continue;

				if(token=="T" || token=="D")
				{
					assert(v.size() == 15 || v.size()==0);
					if(v.size()==0)
						break;

					Vertex v1,v2,v3,n;
					Color surfaceColor, emissionColor;
					v1 = Vertex (v[0],v[1],v[2]);
					v2 = Vertex (v[3],v[4],v[5]);
					v3 = Vertex (v[6],v[7],v[8]);
					surfaceColor = Color (v[9],v[10],v[11], 1); 
					emissionColor = Color(v[12],v[13],v[14],1);

					string tld_token = ( v[12] != 0 ) ? "L" : token;

					Triangle tri = Triangle(v1,v2,v3, tld_token);
					tri.setColor("surface", surfaceColor);
					tri.setColor("emission", emissionColor);

					triangles.push_back(tri);

					v.clear();
				  	continue;
				}
	        	v.push_back(stof(token));
   			}
  		}
  		return triangles;
	}


};