
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


class Ray { 
public:
	Vec Origin;          //ray origin, start point.
	Vec Direction;       //ray direction

	double TMin = 0.0;
	double TMax = MAX_RAY_DISTANCE;

	uint32_t ReflectTimes = 0; //reflect time



	uint32_t depth(){ return ReflectTimes;}

	Ray() = delete;
	Ray(Vec o_, Vec d_, uint32_t _reflect_times) {Origin = _o; Direction= _d; ReflectTimes = _reflect_times;}

	~Ray(){}


	Ray reflect(Vec p, Vec N){

		// O + (-I) =  2*((-I) dot N)*N;
		//        O = I - 2*(I dot N)*N;   //out direction

		inc_depth();

		Vec I = this->d;
		Vec out_dir = (I -  N*(2.0 *(N.dot(I))));//out direction.

		Ray r1(p, out_dir, depth());
		return r1;
	}

	Ray refract(Vec p, Vec N, double eta){
		
		inc_depth();

		Vec I = this->d;
		double k = 1.0 - eta * eta * (1.0 - ((N.dot(I)) * (N.dot(I))));
		Vec out_dir;
		if (k < 0.0){
			out_dir = Vec();
		}
		else{
			out_dir = ((I*eta) - (N*(eta * (N.dot(I)) + sqrt(k))));
		}

		Ray r1(p, out_dir, depth());
		return r1;
	}
}; 


class Camera{
	Vec center;
	double focus = 10;
	Vec normal;
	double reslution[2] = {0.001, 0.001};
	uint32_t size[2]    = {1024, 768};

	//https://github.com/shashwb/C-Ray-Tracer
};


class PerspectiveCamera{
	Vec position;        
	Vec orientation;     
	double nearDistance  = 0.1; 
	double farDistance   = DBL_MAX;  
	double focalDistance = 0.0;
	double heightAngle   = 0.0;  
};



class Light{
public:
	double intensity = 1.0;
	Vec position;

	// Fields common to all subclasses:
	bool         on;     // Whether light is on
	double       intensity;  // Source intensity (0 to 1)
	Vec          color;      // RGB source color

	static void     initClass();

protected:
	Light();              // Makes this abstract
	virtual ~SoLight();
};




class SoDirectionalLight : public SoLight {

	SO_NODE_HEADER(SoDirectionalLight);

	public:
	// Fields (in addition to those in SoLight):
	SoSFVec3f       direction;  // Illumination direction vector

	// Constructor
	SoDirectionalLight();

	SoEXTENDER public:
		// Creates a light source during rendering
		virtual void    GLRender(SoGLRenderAction *action);

	SoINTERNAL public:
		static void     initClass();

	protected:
	virtual ~SoDirectionalLight();
};







class Frame{
private:
	char *buffer = nullptr;

public:
	uint32_t width = 1024;
	uint32_t height = 728;

	uint32_t component = 3;//default RGB888

	Frame() = delete;
	Frame(uint32_t _w, uint32_t _h){ width =_w; height _h; }
	~Frame(){ if(buffer!=nullptr){delete buffer;}};

	void clear(){ memset(ptr, 0, (width*height*component)*sizeof(char)); }

	void write(){};

	void dump(){};
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


	void calcBaryCoords(const Vec &pnt, Vec &bc) const;
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
// smallpt, a Path Tracer by Kevin Beason, 2008 
// Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt 
//        Remove "-fopenmp" for g++ version < 4.2 
// Usage: time ./smallpt 5000 && xv image.ppm 




//###############################################################################

using namespace std;



#define DEBUG


static uint32_t variable_idx_0 = 0;
static uint32_t variable_idx_1 = 0;

bool b_debug = false;
std::string   dump_file_name ="rt_dump.txt";
std::ofstream log_ofs;

void INIT_LOG()
{
	auto ptr = getenv("DEBUG_RT");

	if(ptr!=nullptr){
		b_debug = true;
	}

	if(b_debug==true){
		log_ofs.open(dump_file_name, std::ofstream::out|std::ofstream::trunc);
		log_ofs.flush();
		log_ofs.close();
	}
}


void RPINT_LOG(std::string Func, uint32_t idx ,std::string var_name, double var)
{
	if (b_debug)
	{
		log_ofs.open(dump_file_name, std::ofstream::out|std::ofstream::app);
		//assert(log_ofs.is_open()==true && "Fail to open rt log file.");

		//log_ofs<<Func<<" \t:"  <<idx<<" \t "<< var_name <<": "<<std::hex<< var <<std::endl<<std::endl;
		log_ofs<<Func<<" \t:"  <<idx<<" \t "<< var_name <<": "<< var <<" \t "<<std::hex<< (*((uint64_t*)&var)) <<std::endl;
		//to froce flush.
		log_ofs.flush();
		log_ofs.close();
	}
}







//https://github.com/shashwb/C-Ray-Tracer
//https://www.cs.utexas.edu/~theshark/courses/cs354/lectures/cs354-raytracing_pseudocode.pdf

//#define DOT(V1, V2)
//#define CROSS(V3, V1, V2)


#define RAY_EPSILON   1e-5
#define EPSILON 0.000000000001
#define CROSS(dest,v1,v2) \
	dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
	dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
	dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#define SUB(dest,v1,v2) \
	dest[0]=v1[0]-v2[0]; \
	dest[1]=v1[1]-v2[1]; \
	dest[2]=v1[2]-v2[2];

//double DOT(Vec v1, Vec v2){
//	return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
//}
//	
//Vec CROSS(Vec v1, Vec v2){
//	Vec res;
//
//	res[0]=v1[1]*v2[2]-v1[2]*v2[1]; 
//	res[1]=v1[2]*v2[0]-v1[0]*v2[2]; 
//	res[2]=v1[0]*v2[1]-v1[1]*v2[0];
//	return res;
//}



const uint32_t MAX_RAY_FELECT_CNT = 5;
const double   MAX_RAY_DISTANCE   = 1e20;

const uint32_t MAX_OBJECT_IN_A_BVH_LEAF = 20;


struct Vec {        
	double x, y, z;                  // position, also color (r,g,b) 
	Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; } 

	Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); } 
	Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); } 
	Vec operator*(double b)     const { return Vec(x*b,y*b,z*b); } 
	bool operator==(const Vec &b) const { return ((x==b.x) && (y==b.y) && (z==b.z)); } 

	Vec mult(const Vec &b)      const { return Vec(x*b.x,y*b.y,z*b.z); } 

	double dot(const Vec &b)    const { return x*b.x+y*b.y+z*b.z; }

	//cross 
	Vec operator%(Vec&b)              { return Vec(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x );}
	Vec cross(const Vec &b)     const { return Vec(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x); }

	Vec& norm(){ return *this = *this * (1/sqrt(x*x+y*y+z*z)); } 
	
	double length(){ return sqrt(x*x+y*y+z*z);}
}; 

//see openGL spec.
//genType reflect(genType I, genType N){
//
//	I - 2.0 * dot(N, I) * N
//}
//
//genType refract(genType I, genType N, float eta){
//
//	k = 1.0 - eta * eta * (1.0 - dot(N, I) * dot(N, I));
//	if (k < 0.0)
//		R = genType(0.0);       // or genDType(0.0)
//	else
//		R = eta * I - (eta * dot(N, I) + sqrt(k)) * N;
//}



Vec reflect(Vec I, Vec N){
	//I - 2.0 * dot(N, I) * N;
	return (I -  N*(2.0 *(N.dot(I))));
}

Vec refract(Vec I, Vec N, double eta){
	double k = 1.0 - eta * eta * (1.0 - ((N.dot(I)) * (N.dot(I))));
	if (k < 0.0){
		return Vec();       // or genDType(0.0)
	}
	else{
		return ((I*eta) - (N*(eta * (N.dot(I)) + sqrt(k))));
	}
}





// material types, used in radiance() 
//FIXME: Does a surface will have multi reflect type?
enum Refl_t { 
	DIFF       = 0x0000,
	SPEC       = 0x0001,
	REFR       = 0x0002,


	UNKNOW     = -1
};  


struct Material{
	uint32_t idx;
};

Material materials[100];

//reflecttion, absorb.
class RayStack{
public:
	//
	//See: https://zhuanlan.zhihu.com/p/63177764
	//     https://www.realtimerendering.com/raytracinggems/
	//
	//

	struct Element{
		bool topmost        = false;
		bool odd_parity     = false;
		Material * material = nullptr;
	};

	std::stack<Element> ray_stack;

	void push(Element E){
		ray_stack.push(E);
	};

	void pop(){
		ray_stack.pop();
	};
	// 光线与物体表面相交的问题
	//
	//
	//
	//
	//
	//
};

struct Sphere { 
	double rad;       // radius 
	Vec p;            // position
	Vec e;            // emission, object self light
	Vec c;            // color 
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive) 

	Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_): 
		rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {} 

	Sphere(){};

	double intersect(const Ray &r) const { 
		// returns distance, 0 if nohit 
		// (p'- p)^2 = R^2;
		// ((o+t*d) - p)^2 = R^2;
		// ( t*d + (o-p))^2 - R^2 = 0;
		// d.d*t^2 + 2*d*(o-p)*t + (o-p).(o-p) - R^2 = 0;

		
		// Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 

		Vec op = p-r.o;       //// direction from light_origin to sphere_center
		double b  = (op.dot(r.d));
		double det= b*b-op.dot(op)+rad*rad; 

		double a1 = r.d.dot(r.d);
		double b1 = 2*(r.d.dot(r.o-p));
		double c1 = ((r.o-p).dot(r.o-p)) - rad*rad;
		double delta = b1*b1 - 4*a1*c1;


		if (det<0){
			assert(delta<0);
			return 0; 
		}
		else{
			double eps=1e-4;
			//double r0=sqrt(det); 
			//double t;
			////return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0); 

			//double t0 = b-r0;
			//double t1 = b+r0;

			//if(t0>eps){
			//	return t0;
			//}
			//else if (t1>eps){
			//	return t1;
			//}
			//else{
			//	return 0;
			//}


			double r0 = sqrt(delta);
			
			double t0 = ((0-b1) - r0)/(2*a1);
			double t1 = ((0-b1) + r0)/(2*a1);

			//obviously t0<t1
			if(t0>eps){
				return t0;
			}
			else if (t1>eps){
				return t1;
			}
			else{
				return 0;
			}
		}
	} 
}; 







//object in the 3D space.
Sphere spheres[] = {//Scene: radius, position, emission, color, material 
//#ifdef DEBUG
//	Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left 
//	Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght 
//
//	Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back 
//	//Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(.25,.75,.75),DIFF),//Frnt 
//	Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(), Vec(),            DIFF),//Frnt 
//
//	Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm 
//	Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top 
//
//	//Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr 
//	//Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas 
//	
//	Sphere(15.0, Vec(27,16.5,57),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr 
//
//	Sphere(15.0, Vec(60,16.5,88),       Vec(),Vec(1,1,1)*.999, REFR),//Glas 
//
//	//Sphere(1.0, Vec(60,16.5,96),     Vec(),Vec(0,1,1),  REFR),//Glas 
//
//	Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF) //Lite 
//#else

	Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
	Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
	Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back
	Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(.25,.75,.75),DIFF),//Frnt
	Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm
	Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top

	Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr
	//Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas
	Sphere(16.5,Vec(65,16.5,100),       Vec(),Vec(1,1,1)*.999, REFR),//Glas
	Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF), //Lite
	
	Sphere(5.0, Vec(40,16.5,78),       Vec(),Vec(1,1,1)*.999, DIFF),//a small ball
//#endif

}; 




//struct Line{
//	Vec p[2]; //start end;
//	Vec c[2]; //start end;
//};


#if 0

Vec Cen(50,40.8,-860);
Sphere spheres[] = {//Scene: radius, position, emission, color, material
	// center 50 40.8 62
	// floor 0
	// back  0

	Sphere(1600, Vec(1,0,2)*3000, Vec(1,.9,.8)*1.2e1*1.56*2,Vec(), DIFF), // sun
	Sphere(1560, Vec(1,0,2)*3500,Vec(1,.5,.05)*4.8e1*1.56*2, Vec(),  DIFF), // horizon sun2
	//   Sphere(10000,Cen+Vec(0,0,-200), Vec(0.0627, 0.188, 0.569)*6e-2*8, Vec(.7,.7,1)*.25,  DIFF), // sky
	Sphere(10000,Cen+Vec(0,0,-200), Vec(0.00063842, 0.02001478, 0.28923243)*6e-2*8, Vec(.7,.7,1)*.25,  DIFF), // sky

	Sphere(100000, Vec(50, -100000, 0),  Vec(),Vec(.3,.3,.3),DIFF), // grnd
	Sphere(110000, Vec(50, -110048.5, 0),  Vec(.9,.5,.05)*4,Vec(),DIFF),// horizon brightener
	Sphere(4e4, Vec(50, -4e4-30, -3000),  Vec(),Vec(.2,.2,.2),DIFF),// mountains
	//  Sphere(3.99e4, Vec(50, -3.99e4+20.045, -3000),  Vec(),Vec(.7,.7,.7),DIFF),// mountains snow

	Sphere(26.5,Vec(22,26.5,42),   Vec(),Vec(1,1,1)*.596, SPEC), // white Mirr
	Sphere(13,Vec(75,13,82),   Vec(),Vec(.96,.96,.96)*.96, REFR),// Glas
	Sphere(22,Vec(87,22,24),   Vec(),Vec(.6,.6,.6)*.696, REFR)    // Glas2
};
#endif



//-x,  +x
//+y
//-y


//-z
//
//+z


inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; } 

inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); } 

inline bool intersect(const Ray &r, double &t, int &id){ 
	double n=sizeof(spheres)/sizeof(Sphere);
	double inf=t=1e20; 

	for(int i=int(n);i--;){
		double d=spheres[i].intersect(r);
		if(d!=0 && d<t){
			t=d;
			id=i;
		} 
	}

	return t<inf; 
} 









//return ray_color
Vec radiance(const Ray &r, int depth, unsigned short *Xi){ 

	double t;
	int id;
	if (!intersect(r, t, id)){
		//the light do not hit on anything, just return black.
		return Vec(); // if miss, return black 
	}
	else{
		//the light hit on a object. continue trace the light.

		const Sphere &obj = spheres[id];        //the hit object 
		Vec x=r.o+r.d*t;                        //intersect pos
		Vec n=(x-obj.p).norm();                 //tangent line.

		Vec nl=n.dot(r.d)<0?n:(n*(-1));         //cross or relect.
		Vec f=obj.c; 

		// max refl 
		//max color-component on the oject.
		double p = (f.x>f.y && f.x>f.z) ? f.x : ((f.y>f.z) ? f.y : f.z); 

		if (++depth>MAX_RAY_FELECT_CNT){ 
			//random determin whether continue trace or not.
			//R.R. 
			if (erand48(Xi)<p){ 
				f=f*(1/p);  

				// if a ray have been reflect many-times, 
				// the main component will have more contribution.
				// continue trace the main component "maybe" make more sense.
				//
				// why we need to amplify all component.

				//I think if we "DO NOT" compare random value with max compoent, 
				//and just give it a const threadhold value will both OK.
			}
			else { 
				return obj.e; 
			}

		}


		switch(obj.refl){

			case DIFF:
				{
					// Ideal DIFFUSE reflection
					double r1=2*M_PI*erand48(Xi);
					double r2=erand48(Xi); 
					double r2s=sqrt(r2); 
					Vec w=nl;
					Vec u=((fabs(w.x)>0.1?Vec(0,1,0):Vec(1,0,0))%w).norm(); 
					Vec v= w%u; 
					Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); 
					return obj.e + f.mult(radiance(Ray(x,d),depth,Xi)); 
				}
				break;

			case SPEC:
				{
					// Ideal SPECULAR reflection
					return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi)); 
					//self emission + self_color * other_light_color
				}
				break;

			case REFR:
				{
					Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION 
					bool into = n.dot(nl)>0;                // Ray from outside going in? 
					double nc=1; 
					double nt=1.5; 
					double nnt=into?nc/nt:nt/nc; 
					double ddn=r.d.dot(nl); 
					double cos2t; 
					if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0){    // Total internal reflection 
						return obj.e + f.mult(radiance(reflRay,depth,Xi)); 
					}

					Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm(); 
					double a=nt-nc;
					double b=nt+nc;
					double R0=a*a/(b*b);
					double c = 1-(into?-ddn:tdir.dot(n)); 
					double Re=R0+(1-R0)*c*c*c*c*c;
					double Tr=1-Re;
					double P=.25+.5*Re;
					double RP=Re/P;
					double TP=Tr/(1-P); 

					// Russian roulette 
					return obj.e + f.mult(
							depth>2 ? 
							(
							 erand48(Xi)<P ? 
							 (radiance(reflRay,depth,Xi)*RP)
							 :
							 (radiance(Ray(x,tdir),depth,Xi)*TP)
							) 
							: 
							(radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr)
							); 
				}
				break;

			default:
				assert(0 && "unkonw reflect type.");
		}

#if 0
		if (obj.refl == DIFF){                  // Ideal DIFFUSE reflection 
			double r1=2*M_PI*erand48(Xi), r2=erand48(Xi), r2s=sqrt(r2); 
			Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(), v=w%u; 
			Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); 
			return obj.e + f.mult(radiance(Ray(x,d),depth,Xi)); 
		} else if (obj.refl == SPEC){            // Ideal SPECULAR reflection 
			return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi)); 
		}

		Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION 
		bool into = n.dot(nl)>0;                // Ray from outside going in? 
		double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t; 
		if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0){    // Total internal reflection 
			return obj.e + f.mult(radiance(reflRay,depth,Xi)); 
		}

		Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm(); 
		double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n)); 
		double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P); 
		return obj.e + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette 
					radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) : 
				radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr); 
#endif

	}
} 

double total_t  = 2; //how many seconds we want to render
uint32_t FPS    = 60;//frames per second

void update_scene(uint32_t frame_idx, double delta_t, Vec init_pos, Vec center){

	double r = (2*M_PI)/(1*FPS);  //angular velocity, run a round need 3s
	                             //per frame 60 FPS.
	double theta = r*frame_idx;


	double x = cos(theta);
	double z = sin(theta);


	Vec new_pos;            //move radius
	double rad = init_pos.x - center.x;

	//new_pos.x = init_pos.x + x*rad;
	//new_pos.y = init_pos.y + y*rad;
	//new_pos.z = init_pos.z;

	spheres[9].p.x = center.x + x*rad;
	spheres[9].p.y = center.y;
	spheres[9].p.z = center.z + z*rad;
}

//the main function to render a frame
void render_frame(
	Vec cam_pos /*camera pos*/,
	Vec cam_dir /*camera direction*/,

	Vec *c     /*buffer to save frame*/, 
	const uint32_t w /*image width*/, 
	const uint32_t h /*image hight*/,
	const uint32_t samps /*samples per subpixel*/
){

	//Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir 
	Ray cam(cam_pos, cam_dir.norm()); // cam pos, dir 



	//what is it?
	//Vec cx= Vec(w*.5135/h);
	//Vec cy= (cx%cam.d).norm()*.5135;


	//FOV? 
	// why it is 0.5135 
	Vec cx = Vec(w*.5135/h);
	Vec cy = (cx%cam.d).norm()*0.5135;

	Vec r; //ray_color on each subpixel

#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP 

	for (int y=0; y<h; y++){                       // Loop over image rows 
		fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1)); 

		unsigned short Xi[3]={0,0,y*y*y};       //random seed?
		//unsigned short Xi[3]={0,0,0};         //random seed?

		for (unsigned short x=0; x<w; x++){   // Loop cols 

			uint32_t i = ((h-y-1)*w+x);                    // on which pixel

			for (int sy=0; sy<2; sy++){                    // 2x2 subpixel rows 
				for (int sx=0; sx<2; sx++){                // 2x2 subpixel cols 

					r=Vec(); //init subpixel value as zero.

					//to cast ray on each subpixel
					for (int s=0; s<samps; s++){ 
						double r1=2*erand48(Xi);
						double dx=((r1<1) ? sqrt(r1)-1: 1-sqrt(2-r1));// dx, near to zero. dx<0 | 0>dx

						double r2=2*erand48(Xi);
						double dy=((r2<1) ? sqrt(r2)-1: 1-sqrt(2-r2)); //

						//we get a random direction, cast ray to it.
						Vec d = (
								cx*( ( (sx+.5 + dx)/2 + x)/w - .5) + 
								cy*( ( (sy+.5 + dy)/2 + y)/h - .5) +
								cam.d
								); 

						//we get the ray's (r,g,b)
						Vec t = radiance(Ray(cam.o+d*140,d.norm()),0,Xi); 

						//add all single_ray*ratio color together, to get the subpixel value.
						r = r + t*(1.0/samps);
					} // Camera rays are pushed ^^^^^ forward to start in interior 

					//add subpixel*ratio value to pixel value.
					c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25; 
				} 
			}
		}
	} 
}


void save_frame(
		const Vec *c             /*buffer of a frame*/, 
		const uint32_t w         /*image width*/, 
		const uint32_t h         /*image hight*/,
		const uint32_t frame_idx /*frame idx*/
){


	//we want file name to be:
	// 0000.ppm
	// 0001.ppm
	// 0002.ppm
	// ...
	// 0009.ppm
	// 0010.ppm
	// ...
	// 9999.ppm

	assert((frame_idx <1e4) && "error: we want frame_idx less than 10,000.");

	std::string t_name = std::to_string(frame_idx);

	std::string pre_zero = "0";
	while(t_name.size()<4){
		t_name = "0"+t_name;
	}

#ifdef DEBUG
	std::string file_name = t_name + ".ppm";
#else
	std::string file_name = "./out/" + t_name + ".ppm";
#endif


	FILE *f = fopen(file_name.data(), "w");         // Write image to PPM file. 

	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); //file format
	for (int i=0; i<w*h; i++) {
		fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z)); 
	}

	fflush(f);
	fclose(f);
}

void convert_ppm_to_gif(){
	//convert -delay 0 *ppm movie1.gif

	cout<<"run finish"<<endl;
}


int main_0(int argc, char *argv[]){ 

	uint32_t w=1024; 
	uint32_t h=768;
	uint32_t samps = ((argc==2) ? atoi(argv[1])/4 : 1); // # samples 


	double delta_t = ((double)(1.0))/FPS;

	uint32_t frame_cnt = total_t/delta_t; //how many frame we need to render.

	
	Vec *c=new Vec[w*h]; //buffer to save finnal image. the init value is (0,0,0)
	
	Vec init_pos = spheres[9].p;
	Vec center   = spheres[7].p;

#ifdef DEBUG
	uint32_t frame_idx = 0;
#else
	for(uint32_t frame_idx = 0; frame_idx<frame_cnt; frame_idx++)
#endif
	{
		cout<<"frame_idx:"<< frame_idx <<endl;

		double current_t = (((double)frame_idx) + (0.0))*delta_t; //current timestamp
		double mid_t     = (((double)frame_idx) + (0.5))*delta_t; //mid  timestamp
		double next_t    = (((double)frame_idx) + (1.0))*delta_t; //next timestamp

		memset(c, 0, sizeof(Vec)*w*h);  //init the buffer as black.

		Vec cam_pos = Vec(50,52,295.6);
		Vec cam_dir = Vec(0,-0.042612,-1).norm();

		//update scene
		//update_scene(frame_idx, delta_t, init_pos, center);

		//create BVH


		//render
		render_frame(cam_pos, cam_dir, c, w, h, samps);

		//save
		save_frame(c, w, h, frame_idx);
	}

	//free buffer.
	delete []c;

	//generage a gif
	convert_ppm_to_gif();


	return 0;

}
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////



class Object;

//to save ray-object Hit info.
struct HitInfo{
	Object *object;         //object ptr;
	double  distance;       // hit distance.

	Vec    position;           // position of hit pos
	Vec    normal;             // normal   of hit pos
	Vec    emission;           // emission of hit pos
	Vec    color;              // emission of hit pos
	Refl_t reflection;      // reflection type of hit pos 


	void * info;      //any other information we can save it here.

	HitInfo(){
		//default value for miss, NOT hit on any object
		this->object     = nullptr;
		this->distance   = MAX_RAY_DISTANCE;
		this->position   = Vec(0, 0, 0);
		this->normal     = Vec(0, 0, 0);
		this->emission   = Vec(0, 0, 0);
		this->color      = Vec(0, 0, 0);
		this->reflection = UNKNOW;
	};

	~HitInfo(){
		//default value for miss, NOT hit on any object
		this->object     = nullptr;
		this->distance   = MAX_RAY_DISTANCE;
		this->position   = Vec(0, 0, 0);
		this->normal     = Vec(0, 0, 0);
		this->emission   = Vec(0, 0, 0);
		this->color      = Vec(0, 0, 0);
		this->reflection = UNKNOW;
	};

	HitInfo (const HitInfo& hit){
		//copy constructor
		this->object	    = hit.object;
		this->distance	    = hit.distance;

		this->position	    = hit.position;
		this->normal	    = hit.normal;
		this->emission	    = hit.emission;
		this->color	        = hit.color;
		this->reflection	= hit.reflection;
	}

	HitInfo& operator=(const HitInfo& hit){
		//assign
		this->object	    = hit.object;
		this->distance	    = hit.distance;

		this->position	    = hit.position;
		this->normal	    = hit.normal;
		this->emission	    = hit.emission;
		this->color	        = hit.color;
		this->reflection	= hit.reflection;

		return *this;
	}

	bool operator<(const HitInfo& hit ){
		return ((this->distance) < (hit.distance));
	}

	bool operator>(const HitInfo& hit ){
		return ((this->distance) > (hit.distance));
	}
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

	//Object()= delete;
	//~Object(){};

	//Object(const Object &) = delete;
	//Object &operator=(const Object &) = delete;
};

//a infinit Plane without boundary
//implicit form
//parametirc form
class Plane : public Object{
public:
	Vec    point;   //a point on the plane;
	Vec    normal;  //normal of plane;

	Vec    color;
	Vec    emission;
	Refl_t reflection;

	//
	// NOTE:
	// Since the infinit plane maybe hitted by any ray,
	// after general BVH builder, we need to store the plane in each BVH leaf_node. 
	//
	double min_x() const { return   DBL_MIN;}
	double max_x() const { return   DBL_MAX;}

	double min_y() const { return   DBL_MIN;}
	double max_y() const { return   DBL_MAX;}

	double min_z() const { return   DBL_MIN;}
	double max_z() const { return   DBL_MAX;}

	bool intersect(const Ray &ray,  HitInfo* hit){
		// n dot (p -p0) = 0;
		// n dot (p' - point) = 0;
		// n dot ( (o + t*d) - point) = 0;
		// n dot (o-p) + t* n dot d = 0;
		// t = (n dot (o -p)) / (n dot d);


		double delta = ray.d.dot(this->normal);
		double eps = 1e-4;
		if(delta<eps){
			//the ray is parallel to the plane
			return false;
		}
		else{

			Vec op = (ray.o) - (this->point);
			double upper = normal.dot(op);

			double t = upper/delta;

			if(t >= ray.TMax){
				//too far to hit, just return miss
				return false;
			}
			else{

				hit->object     =  const_cast<Plane *>(this);
				hit->distance   = t;

				hit->position   = ray.o + ray.d*t;
				hit->normal     = this->normal;
				hit->color      = this->color;
				hit->emission   = this->emission;
				hit->reflection = this->reflection;

				return true;
			}
		}
	} 

	//n dot (p -p0) = 0;

	//bool HitTest(const Ray& ray, HitTestResult* result)
	//{
	//	auto denom = Normal * ray.Direction;
	//	if (denom > 0) return false;
	//	auto d = (Normal * ray.Position + Offset) / (-denom);
	//	if (d >= ray.MaxDistance) return false;
	//	result->Shape = const_cast<Plane*>(this);
	//	result->Normal = Normal;
	//	result->Distance = d;
	//	return true;
	//}
};


class Triangle : public Object{
public:
	Vec position[3];        // 3 position. 
	Vec color[3];           // 3 color.
	Vec normal[3];          // 3 normal.

	Vec emission;           // emission, object self light
	Refl_t reflection;      // reflection type (DIFFuse, SPECular, REFRactive)
	
	double min_x() const { return    ((position[0].x < position[1].x) ? ((position[0].x < position[2].x)? position[0].x : position[2].x) :(( position[1].x < position[2].x )? position[1].x : position[2].x));}
	double max_x() const { return    ((position[0].x > position[1].x) ? ((position[0].x > position[2].x)? position[0].x : position[2].x) :(( position[1].x > position[2].x )? position[1].x : position[2].x));}

	double min_y() const { return    ((position[0].y < position[1].y) ? ((position[0].y < position[2].y)? position[0].y : position[2].y) :(( position[1].y < position[2].y )? position[1].y : position[2].y));}
	double max_y() const { return    ((position[0].y > position[1].y) ? ((position[0].y > position[2].y)? position[0].y : position[2].y) :(( position[1].y > position[2].y )? position[1].y : position[2].y));}

	double min_z() const { return    ((position[0].z < position[1].z) ? ((position[0].z < position[2].z)? position[0].z : position[2].z) :(( position[1].z < position[2].z )? position[1].z : position[2].z));}
	double max_z() const { return    ((position[0].z > position[1].z) ? ((position[0].z > position[2].z)? position[0].z : position[2].z) :(( position[1].z > position[2].z )? position[1].z : position[2].z));}


	bool intersect(const Ray &ray,  HitInfo* hit) const {
		//true  for hit
		//false for miss

		//"Moller-Trumbore Alogrithm"

		//to solve:  use Barycentric coordinate system 
		//
		//
		//         p'= (1-u-v)*p0 + u*p1 + v*p2 
		// (o + t*d) = (1-u-v)*p0 + u*p1 + v*p2
		//
		// (o + t*d) = p0 + u*(p1-p0) + v*(p2-p0)
		//
		//  o + t*d  = p0 + u*(p1-p0) + v*(p2-p0)
		//
		// -t*d + u*(p1-p0) + v*(p2-p0) = o-p0
		//
		//                        |t|
		//     [-D, p1-p0, p2-p0] |u| = o - p0;
		//                        |v|
		//     E1 = P1-0;
		//     E2 = P2-0;
		//     T  = O - p0;
		//
		//    Cramer's Rule
		//
		//                        |t|
		//           [-d, E1, E2] |u| = T;
		//                        |v|
		//
		//                        |t|  |T,  E1, E2|
		//                        |u| =|-D, T,  E2|/[-D, E1, E2];
		//                        |v|  |-D, E1, T |
		//
		//
		//
		// this is a "A*x = b" problem
		//
		//     t  = ?
		//     a1 = ?
		//     a2 = ?
		//
		//
		//
		//
		// http://pkuwwt.github.io/scholarship/2014-04-03-ray-triangle-intersection-tests-for-dummies/
		//
		//
		// 
		/*
		*/
# if 0
		Vec O = ray.o;
		Vec D = ray.d;

		Vec E1 = p[1] - p[0];
		Vec E2 = p[2] - p[0];
		Vec T  = O - p[0];

		Vec P = O.cross(E2);
		Vec Q = O.cross(E1);


		double det = P.dot(E1);


		if(abs(det)<EPSILON){
			// the "Denominator" is zero, no means in maths
			//return miss;
			return 0;
		}
		else{

			double t = 0, U = 0, V = 0;

			bool is_frontfacing = true;

			U = T.dot(P);
			V = T.dot(D);

			if(det>EPSILON){
				//front face
				is_frontfacing = true;

				if(U<0 || U>det){
					//return miss;
					return 0;
				}

				if(V<0 || (U+V)>det){
					//return miss;
					return 0;
				}
			}
			else{
				//det< -EPSILON

				//back face
				is_frontfacing = false;
				if(U>0 || U<det){
					//return miss;
					return 0;
				}

				if(V>0 || (U+V)<det){
					//return miss;
					return 0;
				}
			}

			double inv_det = 1.0/det;

			t = Q.dot(E2)*inv_det;
			U = U*inv_det; 
			V = V*inv_det;

			//hit
			return 1;
		}
#endif



		double t, u, v;


		// the precision of the types below can impact speed and accuracy
		// greatly. tweak if you have problems with cracks (or don't).

		double orig[3]   = {ray.o.x, ray.o.y, ray.o.z};
		double dir[3]    = {ray.d.x, ray.d.y, ray.d.z};

		double vert0[3] = {position[0].x, position[0].y, position[0].z};
		double vert1[3] = {position[1].x, position[1].y, position[1].z};
		double vert2[3] = {position[2].x, position[2].y, position[2].z};


		double edge1[3], edge2[3];
		double tvec[3];
		double pvec[3];
		double qvec[3];
		double det;

		// find vectors for two edges sharing vert0
		SUB(edge1, vert1, vert0);
		SUB(edge2, vert2, vert0);

		// begin calculating determinant - also used to calculate U parameter
		CROSS(pvec, dir, edge2);

		// if determinant is near zero, ray lies in plane of triangle
		det = DOT(edge1, pvec);

		if(abs(det)<EPSILON){
			// ray is parallell to the plane of the triangle, return miss;
			return false;  
		}
		else{
			//front facing
			if (det > EPSILON){

				// calculate distance from vert0 to ray origin
				SUB(tvec, orig, vert0);

				// calculate U parameter and test bounds
				u = DOT(tvec, pvec);
				if (u < 0.0 || u > det){
					return false;
				}

				// prepare to test V parameter
				CROSS(qvec, tvec, edge1);

				// calculate V parameter and test bounds
				v = DOT(dir, qvec);
				if (v < 0.0 || u + v > det)
					return false;
			}
			else if(det < -EPSILON){
				//back facing

				// calculate distance from vert0 to ray origin
				SUB(tvec, orig, vert0);

				// calculate U parameter and test bounds
				u = DOT(tvec, pvec);
				if (u > 0.0 || u < det)
					return false;

				// prepare to test V parameter
				CROSS(qvec, tvec, edge1);

				// calculate V parameter and test bounds
				v = DOT(dir, qvec) ;
				if (v > 0.0 || u + v < det)
					return false;
			}

			double inv_det = 1.0 / det;

			// calculate t, ray intersects triangle
			t = DOT(edge2, qvec) * inv_det;
			u = u*inv_det;
			v = v*inv_det;

			double normal[3];
			CROSS(normal,edge1,edge2);

			hit->object     =  const_cast<Triangle *>(this);
			hit->distance   = t;
			hit->position   = ray.o + ray.d*t;
			hit->normal     = Vec(normal[0],normal[1], normal[2]).norm();
			hit->color      = Vec(0.9,0,0);
			hit->emission   = Vec();
			hit->reflection = this->reflection;

			return true;
		}
	}

	Triangle(
			Vec _p0,
			Vec _p1,
			Vec _p2

			//Vec _c0,
			//Vec _c1,
			//Vec _c2
	){
		position[0]=_p0;
		position[1]=_p1;
		position[2]=_p2;

		color[0]=Vec(0.7, 0, 0);
		color[1]=Vec(0.8, 0, 0);
		color[2]=Vec(0.9, 0, 0);

		//c[0]=_c0;
		//c[1]=_c1;
		//c[2]=_c2;

		reflection = DIFF;
	};

	Triangle() = delete;
};

Triangle triangles[]={
	Triangle(Vec(27,16.5,47), Vec(65,16.5,100), Vec(40,16.5,78))
};


void gen_triangle_list(){
	std::vector<Triangle *> list;
	//return 0;
}



class Line : public Object{
public:
	Vec p[2];
	Vec c[2];
};



class Ball: public Object{
public:
	double radious;          // radious
	Vec    position;         // position
	Vec    emission;         // emission, object self light
	Vec    color;            // color 
	Refl_t reflection;       // reflection type (DIFFuse, SPECular, REFRactive) 


	Ball(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_): 
		radious(rad_), position(p_), emission(e_), color(c_), reflection(refl_) {};

	Ball() = delete;



	double min_x() const { return (position.x-radious);}
	double max_x() const { return (position.x+radious);}

	double min_y() const { return (position.y-radious);}
	double max_y() const { return (position.y+radious);}

	double min_z() const { return (position.z-radious);}
	double max_z() const { return (position.z+radious);}

	bool intersect(const Ray &ray,  HitInfo* hit) const {

		// (p'- p)^2 = R^2;
		// ((o+t*d) - p)^2 = R^2;
		// ( t*d + (o-p))^2 - R^2 = 0;
		// d.d*t^2 + 2*d*(o-p)*t + (o-p).(o-p) - R^2 = 0;

		//In math, this is a "a*x^2+ b*y+c=0" problem
		//we can directly use formula to analysis hit or not.

		// t0 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
		// t1 = (-b - sqrt(b^2 - 4*a*c))/(2*a);
		//
		// where:
		// a = d^2;
		// b = 2*d*(o-p);
		// c = (o-p)^2-r^2;

		//The ray direction have be normalized, the d^2 have no effect.
		//thus, we can make some optimization on the formula.


		//Vec op = position - ray.o;
		//double a = 1.0;                          // ray direction have been normalized, this always be 1.0;
		//double b = (ray.d.dot(op));                // 2 have be optimized. and have be pre give a negtive-sign.
		//double c = (op.dot(op)) - radious*radious; //double c = ((r.o-p).dot(r.o-p)) - rad*rad;
		//double delta = b*b - c;                     //this have be optimizd, double delta = (b^2)-(4*a*c);


		Vec op    = (this->position) - ray.o;       //// direction from light_origin to sphere_center
		double b  = (op.dot(ray.d));
		double det= b*b-op.dot(op) + radious*radious ; 

		//double a1 = r.d.dot(r.d);   //equal 1
		//double b1 = 2*(r.d.dot(r.o-position));
		//double c1 = ((r.o-position).dot(r.o-position)) - radious*radious;
		//double delta = b1*b1 - 4*a1*c1;


		if(det<0){
			//miss  
			return false;
		}
		else{
			//hit
			double r0 = sqrt(det);
			
			double t0 = (b - r0);
			double t1 = (b + r0);

			//to find the positive and minimue one.
			//obviously t0<t1;

			double eps=1e-4;
			if(t0>eps){
				//hit on t0;

				hit->object     = const_cast<Ball *>(this);
				hit->distance   = t0;
				hit->position   = ray.o + ray.d*t0;


				Vec n           =  (hit->position - (this->position)).norm();
				Vec nl          =  (n.dot(ray.d)<0)?n:(n*(-1));
				hit->normal     = nl;

				//hit->position = (this->position) + ((hit->normal) * (this->radious-eps));
				//Vec dis0 = (hit->position - (this->position));
				//assert((dis0.length()) >= (this->radious));

				hit->color      = this->color;
				hit->emission   = this->emission;
				hit->reflection = this->reflection;

				return true;
			}
			else if ( t1>eps ){
				//hit on t1;

				hit->object     = const_cast<Ball *>(this);
				hit->distance   = t1;
				hit->position   = ray.o + ray.d*t1;


				Vec n           =  (hit->position -(this->position)).norm();
				Vec nl          =  (n.dot(ray.d)<0)?n:(n*(-1));

				hit->normal     = nl;

				//hit->position = (this->position) + ((hit->normal) * (this->radious-eps));
				//Vec dis0 = (hit->position - (this->position));
				//assert((dis0.length()) >= (this->radious));

				hit->color      = this->color;
				hit->emission   = this->emission;
				hit->reflection = this->reflection;

				return true;
			}
			else {
				//miss
				return false;
			}
		}
	}
	
	
	double intersect(const Ray &r) const { 
		// returns distance, 0 if nohit 
		// (p'- p)^2 = R^2;
		// ((o+t*d) - p)^2 = R^2;
		// ( t*d + (o-p))^2 - R^2 = 0;
		// d.d*t^2 + 2*d*(o-p)*t + (o-p).(o-p) - R^2 = 0;

		
		// Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 

		Vec op    = (this->position) - r.o;       //// direction from light_origin to sphere_center
		double b  = (op.dot(r.d));
		double det= b*b-op.dot(op) + radious*radious ; 

		double a1 = r.d.dot(r.d);   //equal 1
		double b1 = 2*(r.d.dot(r.o-position));
		double c1 = ((r.o-position).dot(r.o-position)) - radious*radious;
		double delta = b1*b1 - 4*a1*c1;


		if (det<0){
			assert(delta<0);
			return 0; 
		}
		else{
			double eps=1e-4;
			double r0=sqrt(det); 
			//double t;
			////return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0); 

			double t0 = b-r0;
			double t1 = b+r0;

			//to find the positive and minimue one.
			//obviously t0<t1;

			if(t0>eps){
				return t0;
			}
			else if (t1>eps){
				return t1;
			}
			else{
				return 0;
			}
		}
	} 
};

//object in the 3D space.
Ball balls[] = {//Scene: radius, position, emission, color, material 
//#ifdef DEBUG
//	Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left 
//	Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght 
//
//	Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back 
//	//Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(.25,.75,.75),DIFF),//Frnt 
//	Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(), Vec(),            DIFF),//Frnt 
//
//	Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm 
//	Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top 
//
//	//Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr 
//	//Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas 
//	
//	Sphere(15.0, Vec(27,16.5,57),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr 
//
//	Sphere(15.0, Vec(60,16.5,88),       Vec(),Vec(1,1,1)*.999, REFR),//Glas 
//
//	//Sphere(1.0, Vec(60,16.5,96),     Vec(),Vec(0,1,1),  REFR),//Glas 
//
//	Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF) //Lite 
//#else

	Ball(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
	Ball(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
	Ball(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back
	Ball(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(.25,.75,.75),DIFF),//Frnt
	Ball(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm
	Ball(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top

	Ball(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, DIFF),//Mirr
	//Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas
	Ball(16.5,Vec(65,16.5,100),       Vec(),Vec(1,1,1)*.999, DIFF),//Glas
	Ball(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF), //Lite
	
	Ball(5.0, Vec(40,16.5,78),       Vec(),Vec(1,1,1)*.999, DIFF),//a small ball
//#endif
}; 

class Ellipsoid: public Object{
	//
	//((x*x)/(a*a))+((y*y)/(b*b))+((z*z)/(c*c))=1
	Vec p;
	Vec r;   //(a,b,c)
};

//Bezier surface
//Bezier Curve and Surface
class BezierSurface: public Object{
public:
	Vec p;
	Vec r;   //(a,b,c)
};

class Curve: public Object{
public:
	Vec p;
};

class Cube: public Object{
	double rad;       // radius 
	Vec p;            // position
	Vec e;            // emission, object self light
	Vec c;            // color 
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive) 
};




#if 0
//aix-align-boundingbox
struct AABB{
	//double min_x;
	//min_x, max_x, min_y, max_y.
	double data[6];
	
	double min_x; double max_x;
	double min_y; double max_y;
	double min_z; double max_z;

	AABB(
	double min_x_ = 0, double max_x_ =0,
	double min_y_ = 0, double max_y_ =0,
	double min_z_ = 0, double max_z_ =0
	):min_x(min_x_),max_x(max_x_),
	min_y(min_y_),max_y(max_y_),
	min_z(min_z_),max_z(max_z_){
		data[0]=min_x_;
		data[1]=max_x_;

		data[2]=min_y_;
		data[3]=max_y_;

		data[4]=min_z_;
		data[5]=max_z_;
	};

	std::vector<Object> obj;

	//to calculate whether a ray hit on the bvh.
	double intersect(Ray r){
		return 0;
	}
};
#endif 


class Scene{
public:
	uint32_t idx;
};



class BVH_NODE{
};


//root
class BVH{
public:
	// AABB : aix-align-boundingbox
	// 0~5: min_x, max_x, min_y, max_y, min_z, max_z
	double aabb[6]={0};
	

	BVH *child0 =nullptr;
	BVH *child1 =nullptr;

	vector<Object *> obj_set;

	BVH(){
		child0 = nullptr;
		child1 = nullptr;
		memset(aabb,0,sizeof(double)*6);

		obj_set.clear();
	};

	~BVH(){

		if(!child0){
			delete child0; 
		}

		if(!child1){
			delete child1; 
		}

		obj_set.clear();
	};

	void init_aabb(double *data){
		memcpy(aabb, data, sizeof(double)*6);
		cout<<"init AABB: \n"
			<< aabb[0] <<"\t"
			<< aabb[1] <<"\t"
			<< aabb[2] <<"\t"
			<< aabb[3] <<"\t"
			<< aabb[4] <<"\t"
			<< aabb[5] <<endl;
	}

	//__attribute((always_inline)) 
	inline bool is_hit_aabb(const Ray& ray){

		//to do aabb intersection analysis.

		//enter and exit on x
		double delta_on_x = (ray.d.x - ray.o.x);
		double t0 = ((aabb[0]- ray.o.x)/delta_on_x);
		double t1 = ((aabb[1]- ray.o.x)/delta_on_x);
		double x_enter = (t0<t1)? t0:t1;
		double x_exit  = (t0>t1)? t0:t1;

		//enter and exit on y
		double delta_on_y = (ray.d.y - ray.o.y);
		double t2 = ((aabb[2]- ray.o.y)/delta_on_y);
		double t3 = ((aabb[3]- ray.o.y)/delta_on_y);
		double y_enter = (t2<t3)? t2:t3;
		double y_exit  = (t2>t3)? t2:t3;

		//enter and exit on z
		double delta_on_z = (ray.d.z - ray.o.z);
		double t4 = ((aabb[4]- ray.o.z)/delta_on_z);
		double t5 = ((aabb[5]- ray.o.z)/delta_on_z);
		double z_enter = (t4<t5)? t4:t5;
		double z_exit  = (t4>t5)? t4:t5;


		double t_enter = (x_enter>y_enter)?((x_enter>z_enter)?x_enter:z_enter):((y_enter>z_enter)?y_enter:z_enter);
		double t_exit  = (x_exit<y_exit)?((x_exit<z_exit)?x_exit:z_exit):((y_exit<z_exit)?y_exit:z_exit);

		if (((t_enter < t_exit) && (t_exit > 0))){
			//hit on the aabb
			return true;
		}
		else{
			//miss on the aabb
			return false;
		}
	}
	
	bool intersect(Ray &ray, unsigned short *Xi ,HitInfo* hit){

		if (!is_hit_aabb(ray)){
			//miss on the BVH node, return miss.
			return false;
		}
		else{
			// hit on the AABB node.

			if( (this->child0 != nullptr )){
				// this is a mid_node, for any mid_node it must have two child.

				// continue traversal all child node, to find the nearest hit point.
				HitInfo hit0;  
				HitInfo hit1; 
				bool res0 = child0->intersect(ray, Xi, &hit0);
				bool res1 = child1->intersect(ray, Xi, &hit1);

				if((res0==false) && (res1==false)){
					//miss on the node
					return false;
				}
				else{

					if(hit0 < hit1){
						//hit on hit0;
						*hit = hit0;
					}
					else{
						//hit on hit1;
						*hit = hit1;
					}
					return true;
				}
			}
			else{
				// this is a BVH leaf_node 

				HitInfo hit_obj;

				//to test all objects in the BVH node
				for(uint32_t i=0;i<this->obj_set.size();i++){

					HitInfo hit_x; 

					Object *obj_x = this->obj_set[i];

					if(obj_x->intersect(ray, &hit_x)){

						// hit on a object, 
						// we need find the nearest one.

						if(hit_x < hit_obj){
							hit_obj = hit_x;
						}
					}
				}

				RPINT_LOG("radiance_1", (variable_idx_0++), "hit_obj.poisition.x", hit_obj.position.x);
				RPINT_LOG("radiance_1", (variable_idx_0++), "hit_obj.poisition.y", hit_obj.position.y);
				RPINT_LOG("radiance_1", (variable_idx_0++), "hit_obj.poisition.z", hit_obj.position.z);

				RPINT_LOG("radiance_1", (variable_idx_0++), "hit_obj.normal.x", hit_obj.normal.x);
				RPINT_LOG("radiance_1", (variable_idx_0++), "hit_obj.normal.y", hit_obj.normal.y);
				RPINT_LOG("radiance_1", (variable_idx_0++), "hit_obj.normal.z", hit_obj.normal.z);

				if(hit_obj.object == nullptr){
					//the ray DO NOT hit on any object in the BVH node, return miss
					hit->color      = Vec();
					return false;
				}
				else{
					// the ray hit on a object,
					// continue trace the ray OR return the emission of object.
				
					//the light hit on a object. continue trace the light.

					//const Sphere &obj = spheres[id];          // the hit object 
					//Vec x =ray.o+ray.d*t;                     // intersect pos
					//Vec n =(x-obj.p).norm();                  // normalized normal of hit pos.

					//Vec nl = n.dot(ray.d)<0?n:(n*(-1));    // properly oriented surface normal.

					Vec f  = hit_obj.color;                  // object color, (BRDF modulator) 

					//max refl 
					//max color-component on the oject.
					double p = (f.x>f.y && f.x>f.z) ? f.x : ((f.y>f.z) ? f.y : f.z); 

					//ray.inc_depth();//increase the depth
					if (ray.depth()>MAX_RAY_FELECT_CNT){ 
						//Russian roulette 

						if(erand48(Xi)>p){ 
							//return obj.e; 
							//NOTE: Hit on the object, return object emission.
							//hit->object     = hit_obj.object;
							//hit->distance   = hit_obj.distance;
							//hit->position   = hit_obj.position;
							//hit->normal     = hit_obj.normal;
							//hit->color      = hit_obj.color;
							//hit->emission   = hit_obj.emission;
							//hit->reflection = hit_obj.reflection;

							*hit = hit_obj;

							return true;
						}
						else{
							//FIXME:
							//NOTE :I think the else-block is un-necessary.
							f=f*(1/p);

							// if a ray have been reflect many-times, 
							// the main component will have more contribution.
							// continue trace the main component "maybe" make more sense.
							//
							// why we need to amplify all component.

							//I think if we "DO NOT" compare random value with max compoent, 
							//and just give it a const threadhold value will both OK.
						}
					}

					switch(hit_obj.reflection){

						case DIFF:
							{
								//
								//For diffuse (not shiny) reflection
								//Sample all lights (non-recursive)
								//Send out additional random sample (recursive)
								//




								uint16_t var0 = Xi[0];
								uint16_t var1 = Xi[1];
								uint16_t var2 = Xi[2];

								double seed0 = ((double)var0);
								double seed1 = ((double)var1);
								double seed2 = ((double)var2);

								RPINT_LOG("radiance_1", (variable_idx_0++), "seed0", seed0);
								RPINT_LOG("radiance_1", (variable_idx_0++), "seed1", seed1);
								RPINT_LOG("radiance_1", (variable_idx_0++), "seed2", seed2);

								// Ideal DIFFUSE reflection

								//step 1: get a random ray
								double r1 = 2*M_PI*erand48(Xi); // random angle around.
								double r2 = erand48(Xi);        //
								double r2s= sqrt(r2);           // Get random distance from center 

								//Vec w=nl;                                                        //z direction
								//Vec u=((fabs(w.x)>0.1? Vec(0,1,0): Vec(1,0,0))%w).norm();        //x direction
								//Vec v= w%u;         //cross.                                     //y direction

								//Vec w = hit_obj.normal;                                               //z direction
								//Vec u =((fabs(w.x)>(0.1)? Vec(0,1,0): Vec(1,0,0)).cross(w)).norm();   //u is perpendicular to w
								//Vec v = w.cross(u);    //cross.                                       //v is perpendicular to u, w
								////random Sampling Unit Hemisphere 1-r2;
								//Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); // d is a random reflect direction.

								////Vec w=nl;
								Vec w = hit_obj.normal;                                               //z direction
								Vec u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(); 
								//Vec v= w%u; 
								Vec v = w.cross(u);    //cross.                                       //v is perpendicular to u, w
								Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); 
					

								RPINT_LOG("radiance_1", (variable_idx_0++), "r1", r1);
								RPINT_LOG("radiance_1", (variable_idx_0++), "r2", r2);
								RPINT_LOG("radiance_1", (variable_idx_0++), "r2s", r2s);

								RPINT_LOG("radiance_1", (variable_idx_0++), "w.x", w.x);
								RPINT_LOG("radiance_1", (variable_idx_0++), "w.y", w.y);
								RPINT_LOG("radiance_1", (variable_idx_0++), "w.z", w.z);

								RPINT_LOG("radiance_1", (variable_idx_0++), "u.x", u.x);
								RPINT_LOG("radiance_1", (variable_idx_0++), "u.y", u.y);
								RPINT_LOG("radiance_1", (variable_idx_0++), "u.z", u.z);

								RPINT_LOG("radiance_1", (variable_idx_0++), "w.x", w.x);
								RPINT_LOG("radiance_1", (variable_idx_0++), "w.y", w.y);
								RPINT_LOG("radiance_1", (variable_idx_0++), "w.z", w.z);

								RPINT_LOG("radiance_1", (variable_idx_0++), "v.x", v.x);
								RPINT_LOG("radiance_1", (variable_idx_0++), "v.y", v.y);
								RPINT_LOG("radiance_1", (variable_idx_0++), "v.z", v.z);

								RPINT_LOG("radiance_1", (variable_idx_0++), "d.x", d.x);
								RPINT_LOG("radiance_1", (variable_idx_0++), "d.y", d.y);
								RPINT_LOG("radiance_1", (variable_idx_0++), "d.z", d.z);


								//return obj.e + f.mult(radiance(Ray(x,d),depth,Xi)); 

								//get a random_ray
								//NOTE: the random ray's depth need "+1".
								uint32_t new_depth = ray.depth() + 1; 
								Ray random_ray = Ray(hit_obj.position, d, new_depth);

								//step 2: Sampling Unit Hemisphere
								//recursively
								HitInfo other;      
								bool res0 = this->intersect(random_ray, Xi, &other); 


								hit->color = hit_obj.emission + f.mult(other.color);

								return true;
							}
							break;

						case SPEC:
							{
								assert(0 && "never go here");

								// Ideal SPECULAR reflection
								//step 1. get a reflect ray
								//NOTE: the reflection ray's depth need "+1".
								uint32_t new_depth = ray.depth() + 1; 
								Ray reflect_ray = ray.reflect(hit_obj.position, hit_obj.normal);

								//step 2. continue trace the reflect_ray
								//recursively
								HitInfo other;
								bool res0 = this->intersect(reflect_ray, Xi, &other); 
								
								hit->color = hit_obj.emission + f.mult(other.color);
								return true;
							}
							break;

						case REFR:
							{

								assert(0 && "un-supported.");

#if 0
								// Ideal dielectric REFRACTION
								//step 1. get a refract ray


								//bool into = n.dot(nl)>0;                // Ray from outside going in? 
								
								Vec n =(hit_obj.position - obj.p).norm();        // normalized normal of hit pos.
								bool into = (n.dot(hit_obj.normal)>0);         // Ray from outside going inside?

								double nc = 1;                          //for vacuum 1.0
								double nt = 1.5;                        //for glass 1.5 ~ 1.7
								double nnt= into?(nc/nt):(nt/nc); 

								// D = B - (D dot N)*N;

								//"Snell's law" / "Snell–Descartes law" / "law of refraction"
								//  sin(I)/sin(O) = (n2/n1)

								Vec D = ray.d;
								Vec N = hit_obj.normal;                 //the real normal
								Vec B = D - N*(N.dot(D));

								Vec T; // the refraction direction.

								Vec D1;  // the reflect direction
								Vec D2;  // the reflect direction


								double ddn   = ray.d.dot(hit_obj.normal); 
								double cos2t = (1-nnt*nnt*(1-ddn*ddn));

								//
								// "Fresnel Reflectance Equation:"
								// Percentage of light is reflected (and what refracted) 
								// from a glass surface based on incident angle (ϴa)





								if(into == true){
									//form outside goto inside. 
									//will both reflection and refraction

									double n  = nc/nt;
									//double F0 = ((n-1)^2)/((n+1)^2);

									double a = (nc - nt);
									double b = (nc + nt);
									double cos_theta = 1.0;
									double c = 1- cos_theta;

									double F0 = (a*a)/(b*b);


									//double Fr = F0 + (1-F0)*((1-cos(ϴa))^5);          //Fresnel reflection.
									double Fr = F0 + (1-F0)*(c*c*c*c*c);                //Fresnel reflection.
									double Tr = 1 - Fr;                                 //Fresnel refraction?

									double P = 0.25 + 0.5 * Fr;                         //probality of reflecting.
									double RP = Fr/P;
									double TP = Tr/(1-P);

									//Glass meterail generally both have reflective and refractive
									Ray reflect_ray = ray.reflect(hit_obj.position, hit_obj.normal);
									Ray refract_ray = ray.refract(hit_obj.position, hit_obj.normal, nnt);


									//both reflection and refraction
									//recursively
									HitInfo reflect_res;
									HitInfo refract_res;  
									bool res0 = intersect(reflect_ray, Xi, &reflect_res); 
									bool res1 = intersect(refract_ray, Xi, &refract_res);

									HitInfo other ;
									other.color = reflect_res.color*RP + refract_res.color*TP;

									//HitInfo res;
									//res.c = obj.e + f.mult(other.color);

									hit->color = hit_obj.emission + f.mult(other.color);
									
									return true;
								}
								else{
									//from inside goto outsied

									//Glass meterail generally both have reflective and refractive
									Ray reflect_ray = ray.reflect(hit_obj.position, hit_obj.normal);

									double theta_i = 0.0;
									double theta_o = 1.0;
									bool TIR = false;
									//(theta_i > theta_o)? true: false;

									//when incident angle greater than the cirtical_angle, it will TIR
									if (TIR==true){
										// if Total internal reflection (TIR)
										// there will be only reflection
										//recursively

										HitInfo reflect_res; 
										bool res0 = intersect(reflect_ray, Xi, &reflect_res); 

										//HitInfo res;
										//res.c = obj.e + f.mult(reflect_res.c);
										
										hit->color = hit_obj.emission + f.mult(reflect_res.color);

										return true;
									}
									else{

										//both reflection and refraction

										//Glass meterail generally both have reflective and refractive
										Ray reflect_ray = ray.reflect(hit_obj.position, hit_obj.normal);
										Ray refract_ray = ray.refract(hit_obj.position, hit_obj.normal, nnt);

										double RP = 0.0, TP = 0.0;

										//both reflection and refraction
										//recursively
										HitInfo reflect_res;
										HitInfo refract_res;  
										bool res0 = intersect(reflect_ray, Xi, &reflect_res); 
										bool res1 = intersect(refract_ray, Xi, &refract_res);

										HitInfo other;
										other.color = reflect_res.color*RP + refract_res.color*TP;

										//HitInfo res;
										//res.c = obj.e + f.mult(other.color);
										
										hit->color = hit_obj.emission + f.mult(other.color);

										return true;
									}
								}
#endif
							}
							break;

						default:
							assert(0 && "unkonw reflect type.");
					}
				}
			}
		}
	}
};


void split_object(
        uint32_t axis, 
        double mid_value, 
      std::vector<Object *> &src_set, 
	  std::vector<Object *> &dst0_set,
	  std::vector<Object *> &dst1_set){

	size_t  cnt = src_set.size();
	assert( cnt!=0 && "invalid src_set");

	size_t half = cnt/2;


	switch(axis){
		case 0://split on x
			{


				for(size_t i=0;i<cnt;i++){

					Object *obj = src_set[i];
					if(i < half){
						dst0_set.push_back(obj);
					}
					else{
						dst1_set.push_back(obj);
					}
				}

			}
			break;

		case 1://split on y
			{

				for(size_t i=0;i<src_set.size();i++){
					Object *obj = src_set[i];
					if((obj->min_x()) >  mid_value){
						dst0_set.push_back(obj);
					}
					else{
						dst1_set.push_back(obj);
					}
				}
			}
			break;

		case 2://split on z
			{
				for(size_t i=0;i<src_set.size();i++){
					Object *obj = src_set[i];
					if((obj->min_x()) >  mid_value){
						dst0_set.push_back(obj);
					}
					else{
						dst1_set.push_back(obj);
					}
				}
			}
			break;

		default:
			assert(0 && "error: unsupport axis.");
	}
}


//std::vector<Object*> obj_set;


bool compare_on_x(Object * i1, Object * i2){
	return ((i1->max_x()) < (i2->max_x()));
}

bool compare_on_y(Object * i1, Object * i2){
	return ((i1->max_y()) < (i2->max_y()));
}

bool compare_on_z(Object * i1, Object * i2){
	return ((i1->max_z()) < (i2->max_z()));
}


BVH * build_BVH(std::vector<Object *> &src_set){

	// 6 boundary value of all those objects on three axis
	//AABB aabb = new AABB( min_x, max_x, min_y, max_y, min_z, max_z);
	
	double t_aabb[6]={
		DBL_MAX, DBL_MIN, 
		DBL_MAX, DBL_MIN, 
		DBL_MAX, DBL_MIN};


	size_t count = src_set.size();
	assert(count!=0 && "error: invalid count!");

	//to find the global longest axis.
	for(size_t i=0;i<count;i++){

		Object * obj = src_set[i];

		double t[6] = {0};//boundray of current object.

		t[0] = obj->min_x();
		t[1] = obj->max_x();

		t[2] = obj->min_y();
		t[3] = obj->max_y();

		t[4] = obj->min_z();
		t[5] = obj->max_z();

		t_aabb[0] = (t_aabb[0]<t[0])?t_aabb[0]:t[0];//min_x
		t_aabb[1] = (t_aabb[1]>t[1])?t_aabb[1]:t[1];//max_x

		t_aabb[2] = (t_aabb[2]<t[2])?t_aabb[2]:t[2];//min_y
		t_aabb[3] = (t_aabb[3]>t[3])?t_aabb[3]:t[3];//max_y

		t_aabb[4] = (t_aabb[4]<t[4])?t_aabb[4]:t[4];//min_z
		t_aabb[5] = (t_aabb[5]>t[5])?t_aabb[5]:t[5];//max_z
	}


	//AABB aabb = new AABB( min_x, max_x, min_y, max_y, min_z, max_z);

	BVH * bvh_ptr = new BVH();

	//init aabb
	bvh_ptr->init_aabb(t_aabb);


	if(count< MAX_OBJECT_IN_A_BVH_LEAF){
		//less than 10 object, as leaf
		bvh_ptr->obj_set = src_set;
		cout<<"\t BVH_LEAF, object count = "<< count<<endl;
	}
	else{

		cout<<"\t MID BVH_LEAF, object count = "<< count<<endl;

		//greater than 10 continue split.
		double x_len = (t_aabb[1] - t_aabb[0]);
		double y_len = (t_aabb[3] - t_aabb[2]);
		double z_len = (t_aabb[5] - t_aabb[4]);


		uint32_t axis    = -1;
		double mid_value = 0;

		if( (x_len >= y_len) && (x_len >= z_len) ){
			//axis = 0;
			//mid_value = (x_len/2) + aabb[0];

			sort(src_set.begin(), src_set.end(), compare_on_x);
		}

		else if( (y_len >= x_len) && (y_len >= z_len) ){
			//axis = 1;
			//mid_value = (y_len/2) + aabb[2];

			sort(src_set.begin(), src_set.end(), compare_on_y);
		}

		else {
			//axis = 2;
			//mid_value = (y_len/2) + aabb[4];
			sort(src_set.begin(), src_set.end(), compare_on_z);
		}


		//split the object to two subset
		size_t half = count/2;
		std::vector<Object *> set0, set1;
		for(size_t i=0;i<count;i++){
			Object *obj = src_set[i];
			if(i <= half){
				set0.push_back(obj);
			}
			else{
				set1.push_back(obj);
			}
		}

		bvh_ptr->child0 = build_BVH(set0);
		bvh_ptr->child1 = build_BVH(set1);
	}

	return bvh_ptr;
}

bool intersect_0(const Ray &r, double &t, int &id){ 
	double inf=t=1e20; 

	size_t n=sizeof(balls)/sizeof(Ball);

	for(size_t i=0;i<n;i++){
		double d=balls[i].intersect(r);
		

		if(d!=0 && d<t){
			t=d;
			id=i;
			//cout<<"ball : "<< i << " hit on : "<<d<<endl;
		}
		else{
		//	if(d==0){
		//		cout<<"ball : "<< i << " miss"<<endl;
		//	}
		//	if(d!=0 && d>=t){
		//		cout<<"ball : "<< i << " miss on infinit."<<endl;
		//	}
		}
	}

	//cout<<"hit on : "<< id << ", t = "<< t<<endl;
	return t<inf; 
} 


//return ray_color
Vec radiance_0(const Ray &r, int depth, unsigned short *Xi){ 
	double t  = DBL_MAX;                          // distance to intersection 
	int id    = -1;                               // id of intersected object 

	if (!intersect_0(r, t, id)){
		//the light do not hit on anything, just return black.
		return Vec(); // if miss, return black 
	}
	else{
		//the light hit on a object. continue trace the light.

		const Ball &obj = balls[id];        //the hit object 
		Vec x=r.o+r.d*t;                        //intersect pos
		Vec n=(x-obj.position).norm();                 //tangent line.

		Vec nl=n.dot(r.d)<0?n:(n*(-1));         //cross or relect.
		Vec f=obj.color; 

		RPINT_LOG("radiance_0", (variable_idx_0++), "hit_pos.x", x.x);
		RPINT_LOG("radiance_0", (variable_idx_0++), "hit_pos.y", x.y);
		RPINT_LOG("radiance_0", (variable_idx_0++), "hit_pos.z", x.z);

		RPINT_LOG("radiance_0", (variable_idx_0++), "nl.x", nl.x);
		RPINT_LOG("radiance_0", (variable_idx_0++), "nl.y", nl.y);
		RPINT_LOG("radiance_0", (variable_idx_0++), "nl.z", nl.z);


		// max refl 
		//max color-component on the oject.
		double p = (f.x>f.y && f.x>f.z) ? f.x : ((f.y>f.z) ? f.y : f.z); 

		if (++depth>MAX_RAY_FELECT_CNT){ 
			//random determin whether continue trace or not.
			//R.R. 
			if (erand48(Xi)<p){ 
				f=f*(1/p);  

				// if a ray have been reflect many-times, 
				// the main component will have more contribution.
				// continue trace the main component "maybe" make more sense.
				//
				// why we need to amplify all component.

				//I think if we "DO NOT" compare random value with max compoent, 
				//and just give it a const threadhold value will both OK.
			}
			else { 
				return obj.emission; 
			}

		}


		switch(obj.reflection){

			case DIFF:
				{
					
					uint16_t var0 = Xi[0];
					uint16_t var1 = Xi[1];
					uint16_t var2 = Xi[2];

					double seed0 = ((double)var0);
					double seed1 = ((double)var1);
					double seed2 = ((double)var2);

					RPINT_LOG("radiance_0", (variable_idx_0++), "seed0", seed0);
					RPINT_LOG("radiance_0", (variable_idx_0++), "seed1", seed1);
					RPINT_LOG("radiance_0", (variable_idx_0++), "seed2", seed2);

					 // Ideal DIFFUSE reflection
					double r1=2*M_PI*erand48(Xi);
					double r2=erand48(Xi); 
					double r2s=sqrt(r2); 

					Vec w=nl;
					Vec u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(); 
					Vec v= w%u; 
					Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); 
					
					RPINT_LOG("radiance_0", (variable_idx_0++), "r1", r1);
					RPINT_LOG("radiance_0", (variable_idx_0++), "r2", r2);
					RPINT_LOG("radiance_0", (variable_idx_0++), "r2s", r2s);

					RPINT_LOG("radiance_0", (variable_idx_0++), "w.x", w.x);
					RPINT_LOG("radiance_0", (variable_idx_0++), "w.y", w.y);
					RPINT_LOG("radiance_0", (variable_idx_0++), "w.z", w.z);

					RPINT_LOG("radiance_0", (variable_idx_0++), "u.x", u.x);
					RPINT_LOG("radiance_0", (variable_idx_0++), "u.y", u.y);
					RPINT_LOG("radiance_0", (variable_idx_0++), "u.z", u.z);

					RPINT_LOG("radiance_0", (variable_idx_0++), "w.x", w.x);
					RPINT_LOG("radiance_0", (variable_idx_0++), "w.y", w.y);
					RPINT_LOG("radiance_0", (variable_idx_0++), "w.z", w.z);
					
					RPINT_LOG("radiance_0", (variable_idx_0++), "v.x", v.x);
					RPINT_LOG("radiance_0", (variable_idx_0++), "v.y", v.y);
					RPINT_LOG("radiance_0", (variable_idx_0++), "v.z", v.z);

					RPINT_LOG("radiance_0", (variable_idx_0++), "d.x", d.x);
					RPINT_LOG("radiance_0", (variable_idx_0++), "d.y", d.y);
					RPINT_LOG("radiance_0", (variable_idx_0++), "d.z", d.z);

					return obj.emission + f.mult(radiance_0(Ray(x,d),depth,Xi)); 

				}
				break;

			case SPEC:
				{
					assert(0 && "never go here.");

					// Ideal SPECULAR reflection
					return obj.emission + f.mult(radiance_0(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi)); 
					//self emission + self_color * other_light_color
				}
				break;

			case REFR:
				{
					assert(0 && "un-supported!");
#if 0
					Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION 
					bool into = n.dot(nl)>0;                // Ray from outside going in? 
					double nc=1; 
					double nt=1.5; 
					double nnt=into?nc/nt:nt/nc; 
					double ddn=r.d.dot(nl); 
					double cos2t; 
					if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0){    // Total internal reflection 
						return obj.e + f.mult(radiance(reflRay,depth,Xi)); 
					}

					Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm(); 
					double a=nt-nc;
					double b=nt+nc;
					double R0=a*a/(b*b);
					double c = 1-(into?-ddn:tdir.dot(n)); 
					double Re=R0+(1-R0)*c*c*c*c*c;
					double Tr=1-Re;
					double P=.25+.5*Re;
					double RP=Re/P;
					double TP=Tr/(1-P); 
					
					// Russian roulette 
					return obj.e + f.mult(
							depth>2 ? 
							(
							 erand48(Xi)<P ? 
							 (radiance(reflRay,depth,Xi)*RP)
							 :
							 (radiance(Ray(x,tdir),depth,Xi)*TP)
							) 
							: 
							(radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr)
							); 
#endif
				}
				break;

			default:
				assert(0 && "unkonw reflect type.");
		}

#if 0
		if (obj.refl == DIFF){                  // Ideal DIFFUSE reflection 
			double r1=2*M_PI*erand48(Xi), r2=erand48(Xi), r2s=sqrt(r2); 
			Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(), v=w%u; 
			Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); 
			return obj.e + f.mult(radiance(Ray(x,d),depth,Xi)); 
		} else if (obj.refl == SPEC){            // Ideal SPECULAR reflection 
			return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi)); 
		}

		Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION 
		bool into = n.dot(nl)>0;                // Ray from outside going in? 
		double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t; 
		if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0){    // Total internal reflection 
			return obj.e + f.mult(radiance(reflRay,depth,Xi)); 
		}

		Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm(); 
		double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n)); 
		double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P); 
		return obj.e + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette 
					radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) : 
			radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr); 
#endif
	}
} 




class RayTraing{
	uint32_t i=0;
};

//the main function to render a frame
void render_frame_0(
	Vec cam_pos /*camera pos*/,
	Vec cam_dir /*camera direction*/,

	Vec *c     /*buffer to save frame*/, 

	BVH *bvh   /*BVH info*/,
	const uint32_t w /*image width*/, 
	const uint32_t h /*image hight*/,
	const uint32_t samps /*samples per subpixel*/

){

	//Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir 
	Ray cam(cam_pos, cam_dir.norm()); // cam pos, dir 



	//what is it?
	//Vec cx= Vec(w*.5135/h);
	//Vec cy= (cx%cam.d).norm()*.5135;


	//FOV? 
	// why it is 0.5135 
	Vec cx = Vec(w*.5135/h);
	Vec cy = (cx%cam.d).norm()*0.5135;

	Vec r; //ray_color on each subpixel


#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP 

	for (int y=0; y<h; y++){                       // Loop over image rows 
		fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1)); 

		unsigned short Xi[3]={0,0,y*y*y};       //random seed?
		//unsigned short Xi[3]={0,0,0};         //random seed?

		for (unsigned short x=0; x<w; x++){   // Loop cols 

			uint32_t i = ((h-y-1)*w+x);                    // on which pixel

			for (int sy=0; sy<2; sy++){                    // 2x2 subpixel rows 
				for (int sx=0; sx<2; sx++){                // 2x2 subpixel cols 

					r=Vec(); //init subpixel value as zero.

					//to cast ray on each subpixel
					for (int s=0; s<samps; s++){ 
						double r1=2*erand48(Xi);
						double dx=((r1<1) ? sqrt(r1)-1: 1-sqrt(2-r1));// dx, near to zero. dx<0 | 0>dx

						double r2=2*erand48(Xi);
						double dy=((r2<1) ? sqrt(r2)-1: 1-sqrt(2-r2)); //

						//we get a random direction, cast ray to it.
						Vec d = (
								cx*( ( (sx+.5 + dx)/2 + x)/w - .5) + 
								cy*( ( (sy+.5 + dy)/2 + y)/h - .5) +
								cam.d
								); 

						Ray ray = Ray(cam.o+d*140,d.norm(), 0);
						
						unsigned short Yi[3]={Xi[0],Xi[1],Xi[2]};       //random seed?

						//we get the ray's (r,g,b)
						//Vec t0 = radiance_0(ray,0,Xi); 

						//auto ptr = getenv("ENALBE_VALUE_CHECK");
						//if(ptr!=nullptr){

							variable_idx_0 = 0;

							HitInfo hit;
							bool res = bvh->intersect(ray, Yi, &hit);
							Vec t = hit.color;
							//assert(t0 == t);
							//t0  = t;
						//}


						//add all single_ray*ratio color together, to get the subpixel value.
						r = r + t*(1.0/samps);

						//exit(-1);

					} // Camera rays are pushed ^^^^^ forward to start in interior 

					//add subpixel*ratio value to pixel value.
					c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*0.25; 
				} 
			}
		}
	} 
}


int main(int argc, char *argv[]){ 

	INIT_LOG();

	uint32_t w=1024; 
	uint32_t h=768;
	uint32_t samps = ((argc==2) ? atoi(argv[1])/4 : 1); // # samples 


	double delta_t = ((double)(1.0))/FPS;

	uint32_t frame_cnt = total_t/delta_t; //how many frame we need to render.

	
	Vec *c=new Vec[w*h]; //buffer to save finnal image. the init value is (0,0,0)
	
	Vec init_pos = spheres[9].p;
	Vec center   = spheres[7].p;

#ifdef DEBUG
	uint32_t frame_idx = 0;
#else
	for(uint32_t frame_idx = 0; frame_idx<frame_cnt; frame_idx++)
#endif
	{
		cout<<"frame_idx:"<< frame_idx <<endl;

		double current_t = (((double)frame_idx) + (0.0))*delta_t; //current timestamp
		double mid_t     = (((double)frame_idx) + (0.5))*delta_t; //mid  timestamp
		double next_t    = (((double)frame_idx) + (1.0))*delta_t; //next timestamp

		memset(c, 0, sizeof(Vec)*w*h);  //init the buffer as black.

		Vec cam_pos = Vec(50,52,295.6);
		Vec cam_dir = Vec(0,-0.042612,-1).norm();

		//update scene
		//update_scene(frame_idx, delta_t, init_pos, center);

		std::vector<Object *> src_set;
		uint32_t num = sizeof(balls)/sizeof(Ball);
		for(uint32_t idx=0;idx<num;idx++){
			src_set.push_back( &(balls[idx]) );
		}
		src_set.push_back( &(triangles[0]) );


		//build bvh;
		BVH *  bvh = build_BVH(src_set);

		//render
		render_frame_0(cam_pos, cam_dir, c, bvh, w, h, samps);
		delete bvh;

		//save
		save_frame(c, w, h, frame_idx);
	}
	
	//free buffer.
	delete []c;

	//generage a gif
	convert_ppm_to_gif();


	return 0;
}
/////////////////////////////////////////////////////////////




void update_step(){
}


#if 0
#endif
