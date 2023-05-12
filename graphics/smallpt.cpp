// smallpt, a Path Tracer by Kevin Beason, 2008 
// Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt 
//        Remove "-fopenmp" for g++ version < 4.2 
// Usage: time ./smallpt 5000 && xv image.ppm 

#include <math.h>   
#include <stdlib.h> 
#include <stdio.h>  
#include <stdint.h>
#include <assert.h>
#include <iostream>
#include <string>
#include <string.h>
#include <vector>
#include <map>
#include <stack>


using namespace std;

//#define DEBUG






const uint32_t MAX_RAY_FELECT_CNT = 5;

struct Vec {        
	double x, y, z;                  // position, also color (r,g,b) 
	Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; } 

	Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); } 
	Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); } 
	Vec operator*(double b) const { return Vec(x*b,y*b,z*b); } 

	Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); } 

	Vec& norm(){ return *this = *this * (1/sqrt(x*x+y*y+z*z)); } 
	double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; }
	Vec cross(const Vec &b) const  { return Vec(y*b.z-z*b.y, z*b.x-x*b.y, x*b.y-y*b.x); }

	Vec operator%(Vec&b){return Vec(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x );}//cross 
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

struct Ray { 
	Vec o;          //ray origin, start point.
	Vec d;          //ray direction
	uint32_t reflect_times; //reflect time

	void inc_depth(){ reflect_times++; }
	uint32_t depth(){ return reflect_times;}

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
}; 

// material types, used in radiance() 
enum Refl_t { 
	DIFF       = 0x0000,
	SPEC       = 0x0001,
	REFR       = 0x0002 
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




			//if((t0 >0) && (t1>0)){
			//	return t0<t1? t0:t1;
			//}
			//else if((t0 >0) && (t1<0)){
			//	return t0;
			//}
			//else if((t0 <0) && (t1>0)){
			//	return t1;
			//}
			//else{
			//	return 0;
			//}

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


Sphere *my_spheres = nullptr;
uint32_t cnt = 0;


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


void render_a_frame(
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
		//unsigned short Xi[3]={0,0,0};       //random seed?

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

#ifdef DEBUG
	FILE *f = fopen("image.ppm", "w");         // Write image to PPM file. 
#else

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

	std::string file_name = "./out/"+ t_name + ".ppm";

	FILE *f = fopen(file_name.data(), "w");         // Write image to PPM file. 

#endif

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


int main(int argc, char *argv[]){ 

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
		update_scene(frame_idx, delta_t, init_pos, center);

		//create BVH


		//render
		render_a_frame(cam_pos, cam_dir, c, w, h, samps);

		//save
		save_frame(c, w, h, frame_idx);
	}

	//free buffer.
	delete []c;

	//generage a gif
	convert_ppm_to_gif();


	return 0;

#ifdef DEBUG

#endif

	uint32_t image_count = 1;
	for(uint32_t image_idx =0; image_idx<image_count ;image_idx++)
	{

		Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir 

		//Vec Pos = spheres[8].p;

		//exit(0);

		//what is it?
		//Vec cx= Vec(w*.5135/h);
		//Vec cy= (cx%cam.d).norm()*.5135;


		//FOV? 
		// why it is 0.5135 
		Vec cx= Vec(w*.5135/h);
		Vec cy= (cx%cam.d).norm()*0.5135;


		Vec *c=new Vec[w*h]; //buffer to save finnal image. the init value is (0,0,0)


		Vec r; //ray_color on each subpixel

#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP 

		for (int y=0; y<h; y++){                       // Loop over image rows 
			fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1)); 

			unsigned short Xi[3]={0,0,y*y*y};       //random seed?
			//unsigned short Xi[3]={0,0,0};       //random seed?

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

#ifdef DEBUG
		FILE *f = fopen("image.ppm", "w");         // Write image to PPM file. 
#else
		std::string file_name_prefix = "image";
		//std::string file_name = "./out/"+ std::to_string(frame_idx) + ".ppm";
		std::string file_name = "./out/"+ std::to_string(0) + ".ppm";
		FILE *f = fopen(file_name.data(), "w");         // Write image to PPM file. 
#endif

		fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 
		for (int i=0; i<w*h; i++) {
			fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z)); 
		}

		delete []c;
	}













	//end
}
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


void angle_to_pos(double theta, double &x, double &y){
	//theta 0~2*PI
	//x pos on x
	//y pos on y

	x=cos(theta);
	y=sin(theta);
}


double rate_to_angle(double r, double step){
	return r*step;
}

void add_a_step(uint32_t step){

	double delta_t = 0;
	double t       = delta_t * step;

}

void updata_ball_position(){

	uint32_t FPS = 60;
	double delta_t = 1.0/FPS;
	double total_t = 3;                    // we want render 3s.

	uint32_t frame_cnt = total_t/delta_t; //how many frame we need to render.



	for(uint32_t i =0; i<frame_cnt; i++){

		double current_t = i*delta_t;         //current timestamp
		double mid_t     = (i+(0.5))*delta_t; //mid  timestamp
		double next_t    = (i+1)*delta_t;     //next timestamp

		double r = (2*M_PI)/3;  //angular velocity, run a round need 3s
		double theta = r*i;

		double x = cos(theta);
		double y = sin(theta);

		Vec init_pos;
		Vec center;

		Vec new_pos;            //move radius
		double rad = init_pos.x - center.x;

		new_pos.x = init_pos.x + x*rad;
		new_pos.y = init_pos.y + y*rad;
		new_pos.z = init_pos.z;
	}
}


class Object;

//to save ray-object Hit info.
struct HitInfo{
	Object *obj;      //object ptr;
	
	double  distance; // hit distance.

	Vec p;            // position of hit pos
	Vec n;            // normal   of hit pos
	Vec e;            // emission of hit pos
	Vec c;            // color    of hit pos
	Refl_t refl;      // reflection type of hit pos 

	void * info;      //any other information we can save it here.
	HitInfo(){
		//default for miss.
		obj = nullptr;
		distance = 0;
	};
};



class Object{
public:
	virtual HitInfo intersect(const Ray &r) const = 0;
	//{
	//	HitInfo hit;
	//	return hit;
	//}

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

class Line : public Object{
	Vec p[2];
};

class Triangle : public Object{
public:
	Vec p[3];// 3 position. 
	Vec c[3];// 3 color.
	Vec n[3];// 3 normal.

	Vec e;            // emission, object self light
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
	
	double min_x() const { return    ((p[0].x < p[1].x) ? ((p[0].x < p[2].x)? p[0].x : p[2].x) :(( p[1].x < p[2].x )? p[1].x : p[2].x));}
	double max_x() const { return    ((p[0].x > p[1].x) ? ((p[0].x > p[2].x)? p[0].x : p[2].x) :(( p[1].x > p[2].x )? p[1].x : p[2].x));}

	double min_y() const { return    ((p[0].y < p[1].y) ? ((p[0].y < p[2].y)? p[0].y : p[2].y) :(( p[1].y < p[2].y )? p[1].y : p[2].y));}
	double max_y() const { return    ((p[0].y > p[1].y) ? ((p[0].y > p[2].y)? p[0].y : p[2].y) :(( p[1].y > p[2].y )? p[1].y : p[2].y));}

	double min_z() const { return    ((p[0].z < p[1].z) ? ((p[0].z < p[2].z)? p[0].z : p[2].z) :(( p[1].z < p[2].z )? p[1].z : p[2].z));}
	double max_z() const { return    ((p[0].z > p[1].z) ? ((p[0].z > p[2].z)? p[0].z : p[2].z) :(( p[1].z > p[2].z )? p[1].z : p[2].z));}


	HitInfo intersect(const Ray &r) const {
		
		//to solve:  use Barycentric coordinate system 
		//
		//
		//         p'= (1-a1-a2)*p0 + a1*p1 + a2*p2 
		//  (o+ t*d) = (1-a1-a2)*p0 + a1*p1 + a2*p2
		//
		//
		//      t*d + a1*(p0-p1)+ a2*(p0-p2) = (p0-o);
		//
		//
		// this is a "A*x = b" problem
		//
		//     t  = ?
		//     a1 = ?
		//     a2 = ?
		//
		//
		//      n0 = (p0-p1).cross (p0-p2); //normal of the triangle 
		//      n1 = d.cross(p0-p1);        //normal of (d) and (p0-p1)
		//      n2 = d.cross(p0-p2);        //normal of (d) and (p0-p2)
		//
		//
		//      (t*d + a1*(p0-p1) + a2*(p0-p2)) dot (N) = (p0-o) dot (N)
		//
		//


		//three edge of the triangle.
		Vec edge_0 = p[1] - p[0];
		Vec edge_1 = p[2] - p[1];
		Vec edge_2 = p[0] - p[2];


		//normal of the triangle plane
		Vec N0 = (edge_0.cross(edge_1)).normal;

		double res = (r.dir).dot(N0);

		double eps=1e-4;
		if( abs(res) < eps ){
			//the ray is parallel with the triangle plane
			return miss;
		}
		else{
			// ray hit with the triangle plane
			// continue analysis whether it hit in the triangle


			//normal of ray and triangle edge
			Vec N1 =         (r.d).cross(edge_1);
			Vec N2 =         (r.d).cross(edge_2);

			//      (t*d + a1*(p0-p1) + a2*(p0-p2)) dot (N0) = (p0-o) dot (N0)
			//      (t*d + a1*(p0-p1) + a2*(p0-p2)) dot (N1) = (p0-o) dot (N1)
			//      (t*d + a1*(p0-p1) + a2*(p0-p2)) dot (N2) = (p0-o) dot (N2)


			double t  = ((p[0]-r.o).dot(N0))/(r.d.dot(N0));
			double a1 = ((p[0]-r.o).dot(N2))/((p[0]-p[1]).dot(N2));
			double a2 = ((p[0]-r.o).dot(N1))/((p[0]-p[2]).dot(N1));

			if( 
				((0<t ) && (t<1 )) &&
				((0<a1) && (a1<1)) &&
				((0<a2) && (a2<1)) &&
				(abs(1-(t + a1 + a2)) <eps)
			){
				HitInfo hit;
				hit.obj = const_cast<Triangle *>(this);

				hit.p    = r.o + (r.d)*t;

				//FIXME: we may need to consider this more.
				hit.n    = N0;
				//hit.n    = n[0]*t + n[1]*a1 + n[2]*a2;

				hit.c    = c[0]*t + c[1]*a1 + c[2]*a2;
				hit.e    = c[0]*t + c[1]*a1 + c[2]*a2;

				hit.refl = refl;

				return hit;
			}
			else{
				return miss;
			}
		}

	}
	
};


class Ball: public Object{
	double rad;       // radious
	Vec p;            // position
	Vec e;            // emission, object self light
	Vec c;            // color 
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive) 

	Ball(Vec _p, double _rad, Vec _e, Vec _c, Refl_t _refl):
		 p(_p), rad(_rad), e(_e), c(_c), refl(_refl){}; 




	double min_x() const { return (p.x-rad);}
	double max_x() const { return (p.x+rad);}

	double min_y() const { return (p.y-rad);}
	double max_y() const { return (p.y+rad);}

	double min_z() const { return (p.z-rad);}
	double max_z() const { return (p.z+rad);}


	HitInfo intersect(const Ray &r) const {
		HitInfo hit;
		hit.obj = const_cast<Ball *>(this);
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


		Vec op = p-r.o;

		//double a = 1.0;           // ray direction have been normalized, this always be 1.0;
		double b = (r.d.dot(op));          // 2 have be optimized. and have be pre give a negtive-sign.
		double c = (op.dot(op)) - rad*rad; //double c = ((r.o-p).dot(r.o-p)) - rad*rad;

		double delta = b*b - c;            //this have be optimizd, double delta = (b^2)-(4*a*c);


		if(delta<0){
			//miss  
			hit.obj = nullptr;
			//return miss;
		}
		else{
			//hit
			double r0 = sqrt(delta);
			
			double t0 = (b - r0);
			double t1 = (b + r0);

			//to find the positive and minimue one.
			//obviously t0<t1;

			double eps = 1e-4;
			if(t0>eps){
				//hit on t0;
				hit.obj = const_cast<Ball *>(this);
			}
			else if (t1 > eps) {
				//hit on t1;
				hit.obj = const_cast<Ball *>(this);
			}
			else {
				//miss
				hit.obj = nullptr;
				//return miss;
			}
		}

		return hit;
	}
};

class Ellipsoid: public Object{
	//
	//((x*x)/(a*a))+((y*y)/(b*b))+((z*z)/(c*c))=1
	Vec p;
	Vec r;   //(a,b,c)
};

//Bezier surface
class BezierSurface: public Object{
	Vec p;
	Vec r;   //(a,b,c)
};

struct Cube{
	double rad;       // radius 
	Vec p;            // position
	Vec e;            // emission, object self light
	Vec c;            // color 
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive) 
};




struct Curve{
public:
}; 




std::vector<Object> object_set;


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

HitInfo miss;

struct BVH{

	// AABB : aix-align-boundingbox
	// 0~5: min_x, max_x, min_y, max_y, min_z, max_z
	double aabb[6]={0};
	

	BVH *child0 =nullptr;
	BVH *child1 =nullptr;

	vector<Object *> obj_set;

	BVH(){
		child0 = nullptr;
		child1 = nullptr;
		obj_set.clear();
	};

	~BVH(){
		//delete child0; 
		//delete child1; 
		//obj_set.clear();
	};

	void init_aabb(double *data){
		memcpy(aabb, data, sizeof(double)*6);
	}

	//bool intersect(Ray r){
	//	return false;
	//};
	
	HitInfo intersect(Ray r){

		//to do aabb intersection analysis.

		//enter and exit on x
		double delta_on_x = (r.d.x - r.o.x);
		double t0 = ((aabb[0]- r.o.x)/delta_on_x);
		double t1 = ((aabb[1]- r.o.x)/delta_on_x);
		double x_enter = (t0<t1)? t0:t1;
		double x_exit  = (t0>t1)? t0:t1;
		
		//enter and exit on y
		double delta_on_y = (r.d.y - r.o.y);
		double t2 = ((aabb[2]- r.o.y)/delta_on_y);
		double t3 = ((aabb[3]- r.o.y)/delta_on_y);
		double y_enter = (t2<t3)? t2:t3;
		double y_exit  = (t2>t3)? t2:t3;

		//enter and exit on z
		double delta_on_z = (r.d.z - r.o.z);
		double t4 = ((aabb[4]- r.o.z)/delta_on_z);
		double t5 = ((aabb[5]- r.o.z)/delta_on_z);
		double z_enter = (t4<t5)? t4:t5;
		double z_exit  = (t4>t5)? t4:t5;


		double t_enter = (x_enter>y_enter)?((x_enter>z_enter)?x_enter:z_enter):((y_enter>z_enter)?y_enter:z_enter);
		double t_exit  = (x_exit<y_exit)?((x_exit<z_exit)?x_exit:z_exit):((y_exit<z_exit)?y_exit:z_exit);

		if((t_enter < t_exit) && (t_exit > 0)){
			//hit on the AABB node.

			if( (this->child0 == nullptr ) && ( this->child1 == nullptr) ){
				// this is a leaf_node 
				// to traversal all object.

				double inf = 0;
				HitInfo hit_obj;

				uint32_t id = -1;
				double    t = 0;

				//to test all objects in the BVH node
				for(uint32_t i=0;i<this->obj_set.size();i++){

					Object *obj = this->obj_set[i];

					HitInfo hit_x = obj->intersect(r);

					if( hit_x.obj != nullptr){
						// hit on a object, 
						// we need find the nearest one.
						if(hit_x.distance < hit_obj.distance){
							//FIXME: to return value data.
							hit_obj = hit_x;
						}
					}
				}

				if(hit_obj.obj == nullptr){
					//return Vec(); //miss, return black 
					//the ray DO NOT hit on any object in the BVH node.
					return miss;
				}
				else{
					// the ray hit on a object,
					// continue trace the ray OR return the color.
#if 1	
					unsigned short *Xi; 
					uint32_t id = 0;
					double t;
					//the light hit on a object. continue trace the light.

					const Sphere &obj = spheres[id];       //the hit object 
					Vec x =r.o+r.d*t;                      //intersect pos
					Vec n =(x-obj.p).norm();               //normalized normal of hit pos.

					Vec nl = n.dot(r.d)<0?n:(n*(-1));      // properly oriented surface normal.
					Vec f  = obj.c;                        // object color, (BRDF modulator) 

					//max refl 
					//max color-component on the oject.
					double p = (f.x>f.y && f.x>f.z) ? f.x : ((f.y>f.z) ? f.y : f.z); 

					r.inc_depth();//increase the depth
					if (r.depth()>MAX_RAY_FELECT_CNT){ 

						if(erand48(Xi)>p){ 
							//return obj.e; 
							//FIXME: Hit on the object, return object emission.
							return hit_obj;
						}
						else{
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

					switch(hit_obj.refl){

						case DIFF:
							{
								//
								//For diffuse (not shiny) reflection
								//Sample all lights (non-recursive)
								//Send out additional random sample (recursive)


								// Ideal DIFFUSE reflection

								//step 1: get a random ray
								double r1 = 2*M_PI*erand48(Xi); // random angle around.
								double r2 = erand48(Xi);        //
								double r2s= sqrt(r2);           // Get random distance from center 

								//Vec w=nl;                                                        //z direction
								//Vec u=((fabs(w.x)>0.1? Vec(0,1,0): Vec(1,0,0))%w).norm();        //x direction
								//Vec v= w%u;         //cross.                                     //y direction

								Vec w = hit_obj.n;                                                    //z direction
								Vec u =((fabs(w.x)>(0.1)? Vec(0,1,0): Vec(1,0,0)).cross(w)).norm();   //u is perpendicular to w
								Vec v = w.cross(u);    //cross.                                       //v is perpendicular to u, w


								//random Sampling Unit Hemisphere 1-r2;
								Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); // d is a random reflect direction.

								//return obj.e + f.mult(radiance(Ray(x,d),depth,Xi)); 

								//get a random_ray
								Ray random_ray = Ray(hit_obj.p, d, r.depth());

								//step 2: Sampling Unit Hemisphere
								HitInfo other      = intersect(random_ray);

								HitInfo res;
								res.c = hit_obj.e + f.mult(other.c);
								return res;
							}
							break;

						case SPEC:
							{
								// Ideal SPECULAR reflection

								//
								//return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi)); 
								//self emission + self_color * other_light_color
								//

								//step 1. get a reflect ray
								Ray reflect_ray = r.reflect(hit_obj.p, hit_obj.n);

								//step 2. continue trace the reflect_ray
								HitInfo other = intersect(reflect_ray);
								
								HitInfo res;
								res.c = obj.e + f.mult(other.c);
								return res;
							}
							break;

						case REFR:
							{
								// Ideal dielectric REFRACTION
								//step 1. get a refract ray

								//FIXME: should we base on different object type?
								//switch(object_type){
								//}

								assert(0 && "un-supported.");


								//bool into = n.dot(nl)>0;                // Ray from outside going in? 
								
								Vec n =(hit_obj.p - obj.p).norm();        // normalized normal of hit pos.
								bool into = (n.dot(hit_obj.n)>0);         // Ray from outside going inside?

								double nc = 1;                          //for vacuum 1.0
								double nt = 1.5;                        //for glass 1.5 ~ 1.7
								double nnt= into?(nc/nt):(nt/nc); 

								// D = B - (D dot N)*N;

								//"Snell's law" / "Snell–Descartes law" / "law of refraction"
								//  sin(I)/sin(O) = (n2/n1)

								Vec D = r.d;
								Vec N = hit_obj.n;                 //the real normal
								Vec B = D - N*(N.dot(D));

								Vec T; // the refraction direction.

								Vec D1;  // the reflect direction
								Vec D2;  // the reflect direction


								double ddn   = r.d.dot(hit_obj.n); 
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
									Ray reflect_ray = r.reflect(hit_obj.p, hit_obj.n);
									Ray refract_ray = r.refract(hit_obj.p, hit_obj.n, nnt);


									//both reflection and refraction
									HitInfo reflect_res = intersect(reflect_ray); 
									HitInfo refract_res = intersect(refract_ray);

									HitInfo other ;
									other.c = reflect_res.c*RP + refract_res.c*TP;

									HitInfo res;
									res.c = obj.e + f.mult(other.c);
									return res;
								}
								else{
									//from inside goto outsied


									//Glass meterail generally both have reflective and refractive
									Ray reflect_ray = r.reflect(hit_obj.p, hit_obj.n);

									double theta_i = 0.0;
									double theta_o = 1.0;
									bool TIR = false;
									//(theta_i > theta_o)? true: false;

									//when incident angle greater than the cirtical_angle, it will TIR
									if (TIR==true){
										// if Total internal reflection (TIR)
										// there will be only reflection

										HitInfo reflect_res = intersect(reflect_ray); 

										HitInfo res;
										res.c = obj.e + f.mult(reflect_res.c);
										return res;
									}
									else{

										//both reflection and refraction



										//Glass meterail generally both have reflective and refractive
										Ray reflect_ray = r.reflect(hit_obj.p, hit_obj.n);
										Ray refract_ray = r.refract(hit_obj.p, hit_obj.n, nnt);

										double RP = 0.0, TP = 0.0;

										//both reflection and refraction
										HitInfo reflect_res = intersect(reflect_ray); 
										HitInfo refract_res = intersect(refract_ray);

										HitInfo other;
										other.c = reflect_res.c*RP + refract_res.c*TP;

										HitInfo res;
										res.c = obj.e + f.mult(other.c);
										return res;
									}
								}


#if 0
								if (cos2t<0){
									//when interface greater than the cirtical angle. it will TIR
									// if Total internal reflection (TIR)
									// same a reflection.
									//return obj.e + f.mult(radiance(reflRay,depth,Xi)); 

									//FIXME: todo
									Hit other = intersect(reflect_ray);
									

									Hit info;
									info.c = (hit_obj.e) + (hit_obj.c).mult(other.c);
									return info;
								}
								else{

									Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm(); 
									double a=nt-nc;
									double b=nt+nc;
									double c = 1-(into?(-ddn):tdir.dot(n)); 


									double R0 = a*a/(b*b);
									double Re = R0+(1-R0)*c*c*c*c*c;

									double Tr=1-Re;
									double P=.25+.5*Re;
									double RP=Re/P;
									double TP=Tr/(1-P); 

									// Russian roulette 
									//FIXME todo
									//		return obj.e + f.mult(
									//				depth>2 ? 
									//				(
									//				 erand48(Xi)<P ? 
									//				 (radiance(reflRay,depth,Xi)*RP)
									//				 :
									//				 (radiance(Ray(x,tdir),depth,Xi)*TP)
									//				) 
									//				: 
									//				(radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr); //need calculate both reflection and refraction
									//				); 

									Hit other;
									if(r.depth()<2){
										//need calculate both reflection and refraction
										//
										//				(radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr)
									}
									else{
										// Russian roulette, random determin continue trace reflection or refraction.
										if(erand48(Xi)<P){
											//				 (radiance(reflRay,depth,Xi)*RP)
										}
										else{
											//				 (radiance(Ray(x,tdir),depth,Xi)*TP)
										}
									}
									

									Hit info;

									info.c = (hit_obj.e) + (hit_obj.c).mult(other.c);
									return info;
								}
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
#endif
				}
			}
			else{
				// this is a mid_node, 
				// continue traversal child node 
				// to find the nearest hit point.

				HitInfo hit0 = child0->intersect(r);
				HitInfo hit1 = child1->intersect(r);

				return (hit0.distance<hit1.distance)?hit0:hit1;
			}
		}
		else{
			//miss on the BVH node, return miss.
			return miss;
		}
	}
};


void split_object(
        uint32_t axis, 
        double mid_value, 
      std::vector<Object *> &src_set, 
	  std::vector<Object *> &dst0_set,
	  std::vector<Object *> &dst1_set){

	switch(axis){
		case 0://split on x
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

		case 1://split on y
			for(size_t i=0;i<src_set.size();i++){
				Object *obj = src_set[i];
				if((obj->min_x()) >  mid_value){
					dst0_set.push_back(obj);
				}
				else{
					dst1_set.push_back(obj);
				}
			}

			break;
		case 2://split on z
			for(size_t i=0;i<src_set.size();i++){
				Object *obj = src_set[i];
				if((obj->min_x()) >  mid_value){
					dst0_set.push_back(obj);
				}
				else{
					dst1_set.push_back(obj);
				}
			}

			break;
		default:
			assert(0 && "error: unsupport axis.");
	}
}


//std::vector<Object*> obj_set;


BVH * build_BVH(std::vector<Object *> src_set){
	//to find the global longest axis.

	double min_x;
	double max_x;

	double min_y;
	double max_y;

	double min_z;
	double max_z;

	
	// 6 boundary value of all those objects on three axis
	double data[6]={};      
	
	// mid_pos value on the tree axis
	double mid_value[3]={}; 

	size_t count = src_set.size();

	std::vector<double> x_low, x_hi;
	std::vector<double> y_low, y_hi;
	std::vector<double> z_low, z_hi;

	double x_low_min, x_low_max;
	double x_hi_min,  x_hi_max;

	double y_low_min, y_low_max;
	double y_hi_min,  y_hi_max;

	double z_low_min, z_low_max;
	double z_hi_min,  z_hi_max;

	double mid_x, mid_y, mid_z;


	//FIXME: to find the mid_value.

	
	for(size_t i=0;i<count;i++){

		Object * obj = src_set[i];

		double t[6] = {0};//boundray of current object.
		
		t[0] = obj->min_x();
		t[1] = obj->max_x();

		t[2] = obj->min_y();
		t[3] = obj->max_y();

		t[4] = obj->min_z();
		t[5] = obj->max_z();


		min_x = (min_x<t[0])?min_x:t[0];
		max_x = (max_x>t[1])?max_x:t[1];

		min_y = (min_y<t[0])?min_y:t[2];
		max_y = (max_y>t[1])?max_y:t[3];

		min_z = (min_z<t[0])?min_z:t[4];
		max_z = (max_z>t[1])?max_z:t[5];
	}


	//AABB aabb = new AABB( min_x, max_x, min_y, max_y, min_z, max_z);

	BVH *bvh_ptr = new BVH;
	//init aabb
	bvh_ptr->init_aabb(data);

	//greater than 10 continue split.
	if(object_set.size()>10){
		
		double x_len = (max_x - min_x);
		double y_len = (max_y - min_y);
		double z_len = (max_z - min_z);

		uint32_t axis = -1;
		double mid_value = 0;
		
		if( (x_len >= y_len) && (x_len >= z_len) ){
			axis = 0;
			mid_value = mid_x;
		}
		if( (y_len >= x_len) && (y_len >= z_len) ){
			axis = 1;
			mid_value = mid_y;
		}
		if( (z_len >= x_len) && (z_len >= y_len) ){
			axis = 2;
			mid_value = mid_x;
		}

		std::vector<Object *> set0, set1;
		split_object(axis, mid_value, src_set, set0, set1);

		bvh_ptr->child0 = build_BVH(set0);
		bvh_ptr->child1 = build_BVH(set1);
	}
	else{
		//leaf
		bvh_ptr->obj_set = src_set;
	}

	return bvh_ptr;
}

BVH *bvh_env=nullptr;

//return ray_color
Vec radiance_0(const Ray &r, int depth, unsigned short *Xi){ 
	double t;                               // distance to intersection 
	int id=0;                               // id of intersected object 

	//Hit hit_info = bvh_env->intersect(r);

	//if(hit_info.obj=nullptr){
	if(true){
		//the light do not hit on anything, just return black.
		return Vec(); // if miss, return black 
	}
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
					Vec u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(); 
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


int main_0(int argc, char *argv[]){ 
	uint32_t w=1024; 
	uint32_t h=768;
	uint32_t samps = ((argc==2) ? atoi(argv[1])/4 : 1); // # samples 

	//for(uint32_t image_idx =0; image_idx<10;image_idx++)
#ifdef DEBUG
	uint32_t image_count = 1;
#else
	uint32_t image_count = 4;
#endif
	for(uint32_t image_idx =0; image_idx<image_count ;image_idx++)
	{

		Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir 

		//Vec Pos = spheres[8].p;

#ifndef DEBUG
		float x = cos(((2*M_PI)/image_count)*image_idx)*18;
		float z = sin(((2*M_PI)/image_count)*image_idx)*18;
		
		spheres[8].p.x =  spheres[8].p.x+x;
		spheres[8].p.z =  spheres[8].p.z+z;
#endif

		//exit(0);

		//what is it?
		//Vec cx= Vec(w*.5135/h);
		//Vec cy= (cx%cam.d).norm()*.5135;


		//FOV? 
		// why it is 0.5135 
		Vec cx= Vec(w*.5135/h);
		Vec cy= (cx%cam.d).norm()*0.5135;

		Vec r; //ray_color on each subpixel

		Vec *c=new Vec[w*h]; //buffer to save finnal image.

		//for(uint32_t i=0;i<(w*h);i++){
		//	c[i].x=1.0;
		//	c[i].y=1.0;
		//	c[i].z=1.0;
		//}


#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP 

		for (int y=0; y<h; y++){                       // Loop over image rows 
			fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1)); 

			unsigned short Xi[3]={0,0,y*y*y};       //random seed?
			//unsigned short Xi[3]={0,0,0};       //random seed?

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


							HitInfo t0 = bvh_env->intersect(Ray(cam.o+d*140,d.norm()));

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

#ifdef DEBUG
		FILE *f = fopen("image.ppm", "w");         // Write image to PPM file. 
#else
		std::string file_name_prefix = "image";
		std::string file_name = "./out/"+ std::to_string(image_idx) + ".ppm";
		FILE *f = fopen(file_name.data(), "w");         // Write image to PPM file. 
#endif

		fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 
		for (int i=0; i<w*h; i++) {
			fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z)); 
		}

		delete []c;
	}
}
/////////////////////////////////////////////////////////////




void update_step(){

}


#if 0
Hit intersect(Ray r, BVH *node){

	if(node->intersect(r)){

		if( (node->child0 == nullptr ) && ( node->child1 == nullptr) ){
			// this is a leaf_node. 
			// to traversal all object.

			double inf = 0;
			Hit hit;

			for(uint32_t i=0;i<node->obj_set.size();i++){
				Object *obj = node->obj_set[i];

				Hit hit_x = obj->intersect(r);

				if(hit_x.distance <hit.distance){
					hit = hit_x;
				}
			}
			return hit;
		}
		else{
			// this is a mid_node, 
			// continue traversal child node.

			Hit hit0 = intersect(r, node->child0);
			Hit hit1 = intersect(r, node->child1);

			return (hit0.distance<hit1.distance)?hit0:hit1;
		}
	}

	//0 for miss
	return hit;
}

#endif
