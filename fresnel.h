#include <iostream>
#include <math.h>
using namespace std;

#ifdef __linux__
        // "Compiled for Linux
#else
        // Windows doesn't define these values by default, Linux does
        #define M_PI 3.141592653589793
        #define INFINITY 1e8
#endif


struct Fresnel{
        float reflection;
        float theta2;
};

float sind(float input){
        return sin(input * M_PI / 180);
}

float cosd(float input){
        return cos(input * M_PI / 180);
}

float tand(float input){
        return tan(input * M_PI / 180);
}

float asind(float input){
        return (asin(input)/ M_PI) *180;
}

float acosd(float input){
        return (acos(input)/ M_PI) *180;
}

Fresnel fresnel(float n1, float n2, float theta1){
        float temp= sind(theta1)*n1/n2;
        Fresnel output;
        if (temp >1){
                output.reflection = 1.0;
                output.theta2 = 90.0;
                return output;
        }
        else {
                output.theta2 = asind(temp);
                float theta2=output.theta2;
                float r1 = pow(((n1* cosd(theta1) - n2*cosd(theta2))/(n1*cosd(theta1)+n2*cosd(theta2))),2);
                float r2 = pow(((n1*cosd(theta2) - n2*cosd(theta1))/(n1*cosd(theta2)+n2*cosd(theta1))),2);
                output.reflection= (r1+r2)/2;
                return output;
        }
}


