
#include <cmath>   
#include <stdlib.h> 
#include <stdio.h>  
#include <stdint.h>
#include <float.h>
#include <limits.h>
#include <assert.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include <map>
#include <stack>


using namespace std;

#define DEBUG


static uint32_t variable_idx_0 = 0;
static uint32_t variable_idx_1 = 0;

bool b_debug = false;
std::string   dump_file_name ="rt_dump.txt";
std::ofstream log_ofs;


class Vec3f;
class Vec2f;
class Vec2s;
class Vec4f;
//class Rotation;
class Matrix;

class Line;
class Plane;
class Cylinder;
class Sphere;


struct Vec {
protected:
	//storage vector value, position(x,y,z) also color (r,g,b) 
	double vec[3] = {0, 0, 0};   


public:
	//Default constructor, all zero
	Vec(){ vec[0] = 0; vec[1] = 0; vec[2] = 0; }

	//Constructor given an array of 3 components
	Vec(double _v[3]){ vec[0] = v[0]; vec[1] = v[1]; vec[2] = v[2]; }

	//Constructor given 3 individual components
	Vec(double x, double y, double z) { vec[0] = x; vec[1] = y; vec[2] = z; }

	//destuctor
	~Vec(){}


	// component-wise scalar multiplication and division operators
	Vec & operator += (const Vec &b) const { vec[0] += b.vec[0]; vec[1] += b.vec[1]; vec[2] += b.vec[2]; return (*this);}
	Vec & operator -= (const Vec &b) const { vec[0] -= b.vec[0]; vec[1] -= b.vec[1]; vec[2] -= b.vec[2]; return (*this);}

	Vec & operator *=(double _v)     const {                       vec[0] *=_v; vec[1] *=_v; vec[2] *=_v; return (*this);}
	Vec & operator /=(double _v)     const { double _d = (1.0/_v); vec[0] *=_d; vec[1] *=_d; vec[2] *=_d; return (*this);}



	friend Vec operator+(const Vec &b) const { return Vec(vec[0]+b.vec[0], vec[1]+b.vec[1], vec[2]+b.vec[2]); } 
	friend Vec operator-(const Vec &b) const { return Vec(vec[0]-b.vec[0], vec[1]-b.vec[1], vec[2]-b.vec[2]); } 

	friend Vec operator*(Vec &a, double _v)   const  { return Vec(a.vec[0]*_v,  a.vec[1]*_v,    a.vec[2]*_v); }
	friend Vec operator*(double _v, Vec &a)   const  { return Vec(a.vec[0]*_v,  a.vec[1]*_v,    a.vec[2]*_v); }
	friend Vec operator/(double _v)   const  { double _d = (1.0/_v); return Vec( vec[0]*_d,  vec[1]*_d,    vec[2]*_d); }

	friend bool operator==(const Vec &b) const { return ((vec[0]==b.vec[0]) && (vec[1]==b.vec[1]) && (vec[2]==b.vec[2])); } 
	friend bool operator!=(const Vec &b) const { return ((vec[0]!=b.vec[0]) || (vec[1]!=b.vec[1]) || (vec[2]!=b.vec[2])); } 



	// Returns dot (inner) product of vector and another vector
	double dot(const Vec &b)    const { return   ((vec[0]*b.vec[0])+(vec[1]*b.vec[1])+(vec[2]*b.vec[2]));}

	//cross 
	//Vec operator%(Vec&b)              { return Vec(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x );}
	// Returns right-handed cross product of vector and another vector
	Vec cross(const Vec &b)     const { return Vec( vec[1]*b.vec[2]-vec[2]*b.vec[1], vec[2]*b.vec[0]-vec[0]*b.vec[2], vec[0]*b.vec[1]-vec[1]*b.vec[0]); }

	// Returns geometric length of vector
	double length const (){ return sqrt( (vec[0]*vec[0]) + (vec[1]*vec[1]) + (vec[2]*vec[2]) );}


	//return the direction of the vector
	Vec norm()      const { double len = length(); return Vec(vec[0]*len, vec[1]*len, vec[2]*len);}

	//Changes vector to be unit length
	Vec& normalize()      { double len = length(); vec[0] *=_len; vec[1] *=_len; vec[2] *=_len; return (*this);}
	
	
	Vec mult(const Vec &b)      const { return Vec(vec[0]*b.vec[0], vec[1]*b.vec[1], vec[2]*b.vec[2]);}
	
	bool    equals(const vec v, double tolerance) const;







	// Returns pointer to array of 3 components
	const double *getValue() const           { return vec; }

	// Returns 3 individual components
	void  getValue(double &x, double &y, double &z) const { x=vec[0]; y=vec[1]; z=vec[2]; };


	// Sets value of vector from array of 3 components
	Vec &   setValue(const float v[3])  { vec[0] = v[0]; vec[1] = v[1]; vec[2] = v[2]; return (*this); }

	// Sets value of vector from 3 individual components
	Vec &   setValue(float x, float y, float z){ vec[0] = x; vec[1] = y; vec[2] = z; return (*this); }
};



class Vec3f{

};



//the abstract for object,
//any detailed object need inherit the Object class
class Object{
public:
	virtual bool intersect(const Ray &ray,  HitInfo* hit) const = 0; 

	virtual double min_x() const = 0;
	virtual double max_x() const = 0;

	virtual double min_y() const = 0;
	virtual double max_y() const = 0;

	virtual double min_z() const = 0;
	virtual double max_z() const = 0;


	void calcBaryCoords(const Vec3f &pnt, Vec &bc) const;
	void useBaryCoords(const double a, const double b, const double c,
			Vec &bc);

	SbColor& getSpecular()         const { return(mat->specular);     }
	bool     isShiny()             const { return(mat->shiny);        }
	bool     isEmissive()          const { return(mat->emissive);     }
	float    getShininess()        const { return(mat->shininess);    }
	float    getTransparency()     const { return(mat->transparency); }
	virtual  float   getEmisAvg()  const { return(mat->emisAvg);      }
	virtual  SbColor getEmission() const { return(mat->emission);     }
	virtual  Vec getNormal()   const {                            };
	Vec  getIntEmission( const Vec &bc) const;
	Vec  getIntNormal(   const Vec &bc) const;
	Vec  getIntDiffColor(const Vec &bc) const;




	void show();
	void show2();



	//Object()= delete;
	//~Object(){};

	//Object(const Object &) = delete;
	//Object &operator=(const Object &) = delete;
};
