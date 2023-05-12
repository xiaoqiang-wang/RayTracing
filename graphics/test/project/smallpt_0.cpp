// smallpt, a Path Tracer by Kevin Beason, 2008 
// Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt 
//        Remove "-fopenmp" for g++ version < 4.2 
// Usage: time ./smallpt 5000 && xv image.ppm 

//#include <math.h>   
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
#include <sstream>
#include <string>
#include <string.h>
#include <vector>
#include <map>
#include <stack>
#include <set>


#include <dirent.h>

//namespace fs = std::filesystem;


using namespace std;

#define DEBUG

//   http://www.4k8k.xyz/article/qjh5606/118117176#_3


static uint32_t variable_idx_0 = 0;
static uint32_t variable_idx_1 = 0;

bool b_debug = false;
std::string   dump_file_name ="rt_dump.txt";
std::ofstream log_ofs;

//  void INIT_LOG()
//  {
//  	auto ptr = getenv("DEBUG_RT");
//  
//  	if(ptr!=nullptr){
//  		b_debug = true;
//  	}
//  
//  	if(b_debug==true){
//  		log_ofs.open(dump_file_name, std::ofstream::out|std::ofstream::trunc);
//  		log_ofs.flush();
//  		log_ofs.close();
//  	}
//  }


//void RPINT_LOG(std::string Func, uint32_t idx ,std::string var_name, double var)
//{
//	if (b_debug)
//	{
//		log_ofs.open(dump_file_name, std::ofstream::out|std::ofstream::app);
//		//assert(log_ofs.is_open()==true && "Fail to open rt log file.");
//
//		//log_ofs<<Func<<" \t:"  <<idx<<" \t "<< var_name <<": "<<std::hex<< var <<std::endl<<std::endl;
//		log_ofs<<Func<<" \t:"  <<idx<<" \t "<< var_name <<": "<< var <<" \t "<<std::hex<< (*((uint64_t*)&var)) <<std::endl;
//		//to froce flush.
//		log_ofs.flush();
//		log_ofs.close();
//	}
//}


#if 1
	#define INIT_LOG()                                                              \
	{                                                                               \
		auto ptr = getenv("DEBUG_RT");                                              \
		if(ptr!=nullptr){                                                           \
			b_debug = true;                                                         \
		}                                                                           \
                                                                                    \
		if(b_debug==true){                                                          \
			log_ofs.open(dump_file_name, std::ofstream::out|std::ofstream::trunc);  \ 
			log_ofs.flush();                                                        \
			log_ofs.close();                                                        \ 
			log_ofs.open(dump_file_name, std::ofstream::out|std::ofstream::app);    \ 
		}                                                                           \
	}



	#define RPINT_LOG(Func, idx ,var_name, var)                                         \
	{                                                                                   \
		if (b_debug)                                                                    \
		{                                                                               \
			assert(log_ofs.is_open()==true && "Fail to open rt log file.");             \
			                                                                            \
			log_ofs<<Func<<" \t:"  <<idx<<" \t "<< var_name                             \
			<<": "<< var <<" \t "<<std::hex<< (*((uint64_t*)&var)) <<std::endl;         \
			log_ofs.flush();                                                            \
		}                                                                               \
	}

#else
	#define INIT_LOG()
	#define RPINT_LOG(Func, idx ,var_name, var)
#endif



//log_ofs.open(dump_file_name, std::ofstream::out|std::ofstream::app);        \


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
const uint32_t MAX_RAY_FELECT_CNT_MUST_RETURN = 10; //a ray reflect 10 times it must return.

const double   MAX_RAY_DISTANCE   = 1e30;

const uint32_t MAX_OBJECT_CNT_IN_A_BVH_LEAF = 4;


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

	double length(){ return sqrt(x*x+y*y+z*z);}

	Vec& norm()     { return *this = (*this) * (1/sqrt(x*x+y*y+z*z)); } 
	Vec& normalize(){ return *this = (*this) * (1/sqrt(x*x+y*y+z*z)); } 
	
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

typedef Vec Vec3f;


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

struct Ray { 
	Vec o;          //ray origin, start point.
	Vec d;          //ray direction


	Vec Origin;          //ray origin, start point.
	Vec Direction;       //ray direction

	double TMax = MAX_RAY_DISTANCE;
	double TMin = 0;

	uint32_t reflect_times; //reflect time

	void inc_depth(){ reflect_times++; }
	uint32_t depth() const { return reflect_times;}

	Ray(){}
	Ray(Vec o_, Vec d_) : o(o_), d(d_),reflect_times(0){} 
	Ray(Vec o_, Vec d_, uint32_t reflect_times_) : o(o_), d(d_),reflect_times(reflect_times_){} 




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


	void dump() const{

		if(b_debug){
			assert(log_ofs.is_open()==true && "Fail to open rt log file.");   

			log_ofs<<"Ray_Info:"<<std::endl
				<<"\t reflect_times : " << reflect_times <<std::endl
				<<"\t o :" << (*((uint64_t *)(&o.x))) <<",\t" << (*((uint64_t *)(&o.y))) <<",\t"<< (*((uint64_t *)(&o.z))) <<std::endl
				<<"\t d :" << (*((uint64_t *)(&d.x))) <<",\t" << (*((uint64_t *)(&d.y))) <<",\t"<< (*((uint64_t *)(&d.z))) <<std::endl
				<<std::endl;
			log_ofs.flush();
		}
	}
}; 


class Camera{
	Vec center;
	double focus = 10;
	Vec normal;
	double reslution[2] = {0.001, 0.001};
	uint32_t size[2]    = {1024, 768};

	//https://github.com/shashwb/C-Ray-Tracer

	enum ViewportMapping {  // What will happen if the camera and viewport aspect ratios differ?

		// The first 3 adjust the viewport to fit the camera.
		CROP_VIEWPORT_FILL_FRAME = 0,  // Draw filled frame around vp.
		CROP_VIEWPORT_LINE_FRAME = 1,  // Draw frame in lines.
		CROP_VIEWPORT_NO_FRAME   = 2,  // Draw no frame.

		ADJUST_CAMERA            = 3,  // Adjust camera to fit viewport.
		LEAVE_ALONE              = 4   // Do nothing. Camera image may become stretched out of proportion
	};






	// Fields
	//SoSFEnum     viewportMapping;// Treatment when aspectRatio not same as viewport's aspectRatio
	Vec          position;       // Location of viewpoint
	//SoSFRotation orientation;    // Orientation (rotation with respect to (0,0,-1) vector)
	double       aspectRatio;    // Ratio of width to height of view plane
	double       nearDistance;   // Distance from viewpoint to view plane
	double       farDistance;    // Distance from viewpoint to far clipping plane
	double       focalDistance;  // Distance from viewpoint to point of focus.
};



class PerspectiveCamera{
	Vec position;        
	Vec orientation;     
	double nearDistance; 
	double farDistance;  
	double focalDistance;
	double heightAngle;  
};


class OrthographicCamera: public Camera{

};



class Light{
public:
	double intensity = 1.0;
	Vec position;
};

//FIXME: Does a surface will have multi reflect type?
enum Refl_t { 
	DIFF       = 0x0000,
	SPEC       = 0x0001,
	REFR       = 0x0002,


	UNKNOW     = -1
};  




class Material{
private:

	uint32_t id;

	enum ReflModel {
		SCHLICK,
		DATA, 
		LIBRARY, 
		BSSRDF, 
		HIERSUBSURF
	};

	//SbColor ambientColor;
	//SbColor specular;
	//SbColor emission;

	double rough1 = 1.0, rough2 = 1.0;  // roughness  (0=mirror, 1=diffuse (default))
	double sym1 = 1.0,   sym2 = 1.0;    // anisotropy, 0=aniisotropic, 1=isotropic (symmetrical, default))

	double R0 = 1.0;                    // reflection at normal incidence; Default: 1, -1 means auto-alculated based on index of refractions


	double transparency; // 0 (opaque) - 1 (clear)
	double ior1;         // OUTSIDE index of refraction Default: 1.0 (air)
	double ior2;         // INSIDE index of refraction Default: 1.5 (glass)
	double atten;        // attenuation (fraction) at 1m when going through inside, Default (0.98, 0.99, 0.98) from specularColor ?


	double emisEcc;    // Light eccentricity (default 1 = diffuse)
	double minDist;  // minimum distance for intersection (default 0.0001)




	double g, s, d;           // internal (private, automatically calculated)

	Vec3f ambientColor;
	Vec3f specular;
	Vec3f emission;

	int shiny;
	int emissive;
	float shininess;
	//float transparency;
	float emisAvg;

public:
	Material(){};
	~Material(){};
};





class Texture{
public:
	uint32_t id;
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
	double t;                               // distance to intersection
	int id=0;                               // id of intersected object
	if (!intersect(r, t, id)) return Vec(); // if miss, return black
	const Sphere &obj = spheres[id];        // the hit object
	Vec x=r.o+r.d*t, n=(x-obj.p).norm(), nl=n.dot(r.d)<0?n:n*-1, f=obj.c;
	double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl
	if (++depth>5) if (erand48(Xi)<p) f=f*(1/p); else return obj.e; //R.R.
	if (obj.refl == DIFF){                  // Ideal DIFFUSE reflection
		double r1=2*M_PI*erand48(Xi), r2=erand48(Xi), r2s=sqrt(r2);
		Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(), v=w%u;
		Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
		return obj.e + f.mult(radiance(Ray(x,d),depth,Xi));
	} else if (obj.refl == SPEC)            // Ideal SPECULAR reflection
		return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi));
	Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION
	bool into = n.dot(nl)>0;                // Ray from outside going in?
	double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;
	if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
		return obj.e + f.mult(radiance(reflRay,depth,Xi));
	Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
	double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
	double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
	return obj.e + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette
				radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
			radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
}



int main_0(int argc, char *argv[]){
	int w=1024, h=768, samps = argc==2 ? atoi(argv[1])/4 : 1; // # samples
	Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir
	Vec cx=Vec(w*.5135/h), cy=(cx%cam.d).norm()*.5135, r, *c=new Vec[w*h];
#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
	for (int y=0; y<h; y++){                       // Loop over image rows
		fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
		for (unsigned short x=0, Xi[3]={0,0,y*y*y}; x<w; x++)   // Loop cols
			for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++)     // 2x2 subpixel rows
				for (int sx=0; sx<2; sx++, r=Vec()){        // 2x2 subpixel cols
					for (int s=0; s<samps; s++){
						double r1=2*erand48(Xi), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
						double r2=2*erand48(Xi), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);
						Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
							cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
						r = r + radiance(Ray(cam.o+d*140,d.norm()),0,Xi)*(1./samps);
					} // Camera rays are pushed ^^^^^ forward to start in interior
					c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
				}
	}
	FILE *f = fopen("image.ppm", "w");         // Write image to PPM file.
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
	for (int i=0; i<w*h; i++)
		fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
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


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////



class Object;
class Line;
class LineSegment;
class Plane;
class Triangle;
class Ball;
class Box;
class BezierSurface;




//to save ray-object Hit info.
struct HitInfo{
	Object *object;         //object ptr;
	double  distance;       // hit distance.
	uint64_t max_reflect_times;//max reflect/refract_times.

	Vec    position;           // position of hit pos
	Vec    normal;             // normal   of hit pos
	Vec    emission;           // emission of hit pos
	Vec    color;              // emission of hit pos
	Refl_t reflection;      // reflection type of hit pos 


	void * info;      //any other information we can save it here.

	HitInfo(){
		//default value for miss, NOT hit on any object
		this->object            = nullptr;
		this->distance          = MAX_RAY_DISTANCE;
		this->max_reflect_times = 0;
		this->reflection        = UNKNOW;
	};

	~HitInfo(){
		this->object            = nullptr;
		this->distance          = MAX_RAY_DISTANCE;
		this->max_reflect_times = 0;
		this->reflection        = UNKNOW;
	};

	HitInfo (const HitInfo& hit){
		//copy constructor
		this->object	       = hit.object;
		this->distance	       = hit.distance;
		this->max_reflect_times = hit.max_reflect_times;

		this->position	       = hit.position;
		this->normal	       = hit.normal;
		this->emission	       = hit.emission;
		this->color	           = hit.color;
		this->reflection	   = hit.reflection;
	}

	HitInfo& operator=(const HitInfo& hit){
		//assign
		this->object	       = hit.object;
		this->distance	       = hit.distance;
		this->max_reflect_times = hit.max_reflect_times;

		this->position	       = hit.position;
		this->normal	       = hit.normal;
		this->emission	       = hit.emission;
		this->color	           = hit.color;
		this->reflection	   = hit.reflection;

		return *this;
	}

	bool operator<(const HitInfo& hit ){
		return ((this->distance) < (hit.distance));
	}

	bool operator>(const HitInfo& hit ){
		return ((this->distance) > (hit.distance));
	}

	void dump() const{
		if(b_debug){
			assert(log_ofs.is_open()==true && "Fail to open rt log file.");   

			log_ofs<<"Hit_Info:"<<"\n"
				<<"\t object : " << object<<"\n"
				<<"\t reflection : " << reflection <<std::endl
				<<"\t max_reflect_times : " << max_reflect_times <<std::endl

				<<"\t distance :" << distance << std::hex << *((uint64_t *)(&distance))<<std::endl
				<<"\t position :" << (*((uint64_t *)(&position.x))) <<",\t" << (*((uint64_t *)(&position.y))) <<",\t"<< (*((uint64_t *)(&position.z))) <<std::endl
				<<"\t normal   :" << (*((uint64_t *)(&normal.x)))   <<",\t" << (*((uint64_t *)(&normal.y)))   <<",\t"<< (*((uint64_t *)(&normal.z)))   <<std::endl
				<<"\t emission :" << (*((uint64_t *)(&emission.x))) <<",\t" << (*((uint64_t *)(&emission.y))) <<",\t"<< (*((uint64_t *)(&emission.z))) <<std::endl
				<<"\t color : "   << (*((uint64_t *)(&color.x)))    <<",\t" << (*((uint64_t *)(&color.y)))    <<",\t"<< (*((uint64_t *)(&color.z)))    <<std::endl
				<<std::endl;

			log_ofs.flush();
		}
	}

};


//the abstract for object,
//any detailed object need inherit the Object class
class Object{
private:
	uint32_t id = 0;

public:
	virtual double intersect(const Ray &ray) const = 0;
	virtual void get_HitInfo(const Ray &ray, const double t, HitInfo *hit) const = 0;



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




//////////////////////////////////////////////////////////////////////////////
//
//  Class: Point
//
//  Represents a small point in 3D space.
//
//////////////////////////////////////////////////////////////////////////////
class Point: public Object{

public:
	double radious;          // radious
	Vec    position;         // position
	Vec    center;           // position
	Vec    emission;         // emission, object self light
	Vec    color;            // color 
	Refl_t reflection;       // reflection type (DIFFuse, SPECular, REFRactive) 


public:

	Point(){};
	Point(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_): 
		radious(rad_), position(p_), emission(e_), color(c_), reflection(refl_) {};

	~Point(){}; 

	void  setCenter(const Vec3f &c) { center = c; }
	void  setRadius(double r)       { radious = r;}
	void  setValue(const Vec3f &c, double r ){ center = c; radious = r;}

	const Vec3f & getCenter() const { return center;}
	double getRadius()        const { return radious;}




	double min_x() const { return (position.x-radious);}
	double max_x() const { return (position.x+radious);}

	double min_y() const { return (position.y-radious);}
	double max_y() const { return (position.y+radious);}

	double min_z() const { return (position.z-radious);}
	double max_z() const { return (position.z+radious);}


	void get_HitInfo(const Ray &ray, const double t, HitInfo *hit) const {

		hit->object            = const_cast<Point *>(this);

		hit->distance          = t;
		hit->max_reflect_times = ray.reflect_times;
		hit->position          = ray.o + ray.d*t;

		Vec n                  =  (hit->position -(this->position)).norm();
		Vec nl                 =  (n.dot(ray.d)<0)?n:(n*(-1));
		hit->normal            = nl;

		hit->color             = this->color;
		hit->emission          = this->emission;
		hit->reflection        = this->reflection;
	}

	double intersect(const Ray &ray) const {

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
			return 0;

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
				return t0;
			}
			else if ( t1>eps ){
				return t1;
			}
			else {
				return 0;
			}
		}
	}
	
	bool  intersect(const Line &l, Vec3f &enter, Vec3f &exit) const;
};

//////////////////////////////////////////////////////////////////////////////
//
//  Class: Line
//
//  Represents a directed line in 3D space.
//
//////////////////////////////////////////////////////////////////////////////
class Line : public Object{
private:
	// Parametric description:
	//  l(t) = pos + t * dir
	Vec3f pos;
	Vec3f dir;
public:
	void setValue(const Vec3f &_p0, const Vec3f &_p1) { pos = _p0; dir = (_p0 - _p1); dir.normalize();}


	Line(const Vec3f &_p0, const Vec3f &_p1, const Vec3f &_c){ setValue(_p0, _p1); }




	// Accessors
	const Vec3f & getPosition()  const { return pos; }
	const Vec3f & getDirection() const { return dir; }


	// Returns the closest point on the line to the given point.
	Vec3f  getClosestPoint(const Vec3f &point) const;


	// Find closest points between the two lines. Return FALSE if they are parallel, otherwise return TRUE.
	bool    getClosestPoints(const Line  &line2, Vec3f &ptOnThis, Vec3f &ptOnLine2 ) const;




	////////////////////////////////////////////////////////////////////////
	//
	// Description:
	//    Intersects the line with the triangle formed bu v0, v1, v2.
	//    Returns TRUE if there is an intersection. If there is an
	//    intersection, barycentric will contain the barycentric
	//    coordinates of the intersection (for v0, v1, v2, respectively)
	//    and front will be set to TRUE if the ray intersected the front
	//    of the triangle (as defined by the right hand rule).
	//
	// Use: internal
	//
	////////////////////////////////////////////////////////////////////////
	bool intersect(const Triangle &triangle, Vec3f &intersection, Vec3f &barycentric, bool &front) const;




	////////////////////////////////////////////////////////////////////////
	//
	// Description:
	//    Intersect the line with a line segment.  The line is augmented with
	//    an angle to form a cone.  The line segment must lie within pickAngle
	//    of the line.  If an intersection occurs, the intersection point is
	//    placed in intersection.
	//
	////////////////////////////////////////////////////////////////////////
	bool intersect(const LineSegment &linesegment,  Vec3f &intersection) const;



	////////////////////////////////////////////////////////////////////////
	//
	// Description:
	//    Intersect the line with a 3D box.  The line is intersected with the
	//    twelve triangles that make up the box.
	//
	//
	////////////////////////////////////////////////////////////////////////
	bool intersect(const Box &box, Vec3f &enter, Vec3f &exit) const;

	
	
	double intersect(const Ray &ray) const
	{
		return 0;
	}

	void get_HitInfo(const Ray &ray, const double t, HitInfo *hit) const {

	}

};

//////////////////////////////////////////////////////////////////////////////
//
//  Class: Linesegment
//
//  Represents a LineSegment in 3D space.
//
//////////////////////////////////////////////////////////////////////////////
class LineSegment : public Object{
private:
	Vec3f start, end;
	Vec3f color;
	double thickness;
public:
	LineSegment() = delete;
	~LineSegment(){};
	LineSegment(const Vec3f &_start, const Vec3f &_end, const Vec3f &_color, double _thickness){

		start     = _start;
		end       = _end;
		color     = _color;
		thickness = _thickness;
	}


	void setPos(const Vec3f &_start, const Vec3f &_end){
		start     = _start;
		end       = _end;
	}

	void setColor(const Vec3f &_color){ color  = _color;}
	void setThickness(const double _t){ thickness = _t;};




	double min_x() const { return  ((start.x< end.x)? (start.x) : (end.x)); }
	double max_x() const { return  ((start.x> end.x)? (start.x) : (end.x)); }

	double min_y() const { return  ((start.y< end.y)? (start.y) : (end.y)); } 
	double max_y() const { return  ((start.y> end.y)? (start.y) : (end.y)); } 

	double min_z() const { return  ((start.z< end.z)? (start.z) : (end.z)); } 
	double max_z() const { return  ((start.z> end.z)? (start.z) : (end.z)); } 



	double intersect(const Ray &ray) const{


		return 0;
	}

	void get_HitInfo(const Ray &ray, const double t, HitInfo *hit) const {

	}
};





//////////////////////////////////////////////////////////////////////////////
//
//  Class: Triangle
//
//  Represents a Triangle in 3D space.
//
//////////////////////////////////////////////////////////////////////////////
class Triangle : public Object {
private:
	Vec position[3];        // 3 position. 
	Vec color[3];           // 3 color.
	Vec normal[3];          // 3 normal.

	Vec normal_p;           // normal of the triangle plane.

	Vec emission;           // emission, object self light
	Refl_t reflection;      // reflection type (DIFFuse, SPECular, REFRactive)


public:
	Triangle(
			Vec _p0,
			Vec _p1,
			Vec _p2,

			Vec _c0,
			Vec _c1,
			Vec _c2,

			Vec _e,
			Refl_t _ref
	){
		position[0] = _p0;
		position[1] = _p1;
		position[2] = _p2;

		//color[0]=Vec(0.7, 0, 0);
		//color[1]=Vec(0.8, 0, 0);
		//color[2]=Vec(0.9, 0, 0);

		color[0]   = _c0;
		color[1]   = _c1;
		color[2]   = _c2;

		//normal_p can be caulated.

		emission   = _e;
		reflection = _ref;
	};

	Triangle() = delete;
	~Triangle(){}


	//////////////
	double min_x() const { return    ((position[0].x < position[1].x) ? ((position[0].x < position[2].x)? position[0].x : position[2].x) :(( position[1].x < position[2].x )? position[1].x : position[2].x));}
	double max_x() const { return    ((position[0].x > position[1].x) ? ((position[0].x > position[2].x)? position[0].x : position[2].x) :(( position[1].x > position[2].x )? position[1].x : position[2].x));}

	double min_y() const { return    ((position[0].y < position[1].y) ? ((position[0].y < position[2].y)? position[0].y : position[2].y) :(( position[1].y < position[2].y )? position[1].y : position[2].y));}
	double max_y() const { return    ((position[0].y > position[1].y) ? ((position[0].y > position[2].y)? position[0].y : position[2].y) :(( position[1].y > position[2].y )? position[1].y : position[2].y));}

	double min_z() const { return    ((position[0].z < position[1].z) ? ((position[0].z < position[2].z)? position[0].z : position[2].z) :(( position[1].z < position[2].z )? position[1].z : position[2].z));}
	double max_z() const { return    ((position[0].z > position[1].z) ? ((position[0].z > position[2].z)? position[0].z : position[2].z) :(( position[1].z > position[2].z )? position[1].z : position[2].z));}


	double intersect(const Ray &ray) const {
		
		double t, u, v;

		// the precision of the types below can impact speed and accuracy
		// greatly. tweak if you have problems with cracks (or don't).

		double orig[3]   = {ray.o.x, ray.o.y, ray.o.z};
		double dir[3]    = {ray.d.x, ray.d.y, ray.d.z};

		double vert0[3] = {position[0].x, position[0].y, position[0].z};
