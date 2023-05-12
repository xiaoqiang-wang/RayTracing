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
			return 0;
		}
		else{
			//front facing
			if (det > EPSILON){

				// calculate distance from vert0 to ray origin
				SUB(tvec, orig, vert0);

				// calculate U parameter and test bounds
				u = DOT(tvec, pvec);
				if (u < 0.0 || u > det){
					return 0;
				}

				// prepare to test V parameter
				CROSS(qvec, tvec, edge1);

				// calculate V parameter and test bounds
				v = DOT(dir, qvec);
				if (v < 0.0 || u + v > det)
					return 0;
			}
			else if(det < -EPSILON){
				//back facing

				// calculate distance from vert0 to ray origin
				SUB(tvec, orig, vert0);

				// calculate U parameter and test bounds
				u = DOT(tvec, pvec);
				if (u > 0.0 || u < det)
					return 0;

				// prepare to test V parameter
				CROSS(qvec, tvec, edge1);

				// calculate V parameter and test bounds
				v = DOT(dir, qvec) ;
				if (v > 0.0 || u + v < det)
					return 0;
			}

			double inv_det = 1.0 / det;

			// calculate t, ray intersects triangle
			t = DOT(edge2, qvec) * inv_det;
			// u = u*inv_det;
			// v = v*inv_det;

			return t;
		}
	}

	void get_HitInfo(const Ray &ray, const double t, HitInfo *hit) const {

		double orig[3]   = {ray.o.x, ray.o.y, ray.o.z};
		double dir[3]    = {ray.d.x, ray.d.y, ray.d.z};

		double vert0[3] = {position[0].x, position[0].y, position[0].z};
		double vert1[3] = {position[1].x, position[1].y, position[1].z};
		double vert2[3] = {position[2].x, position[2].y, position[2].z};


		double edge1[3], edge2[3];
		
		// find vectors for two edges sharing vert0
		SUB(edge1, vert1, vert0);
		SUB(edge2, vert2, vert0);

		double normal[3];
		CROSS(normal,edge1,edge2);


		hit->object     =  const_cast<Triangle *>(this);
		hit->distance   = t;
		hit->position   = ray.o + ray.d*t;
		hit->normal     = Vec(normal[0],normal[1], normal[2]).norm();
		//hit->color      = Vec(1,1,1)*.999;
		hit->color      = Vec(0.9,0,0);
		hit->emission   = Vec();
		hit->reflection = this->reflection;
	}

};

Triangle triangles[]={
	Triangle(Vec(27,16.5,47), Vec(65,16.5,100), Vec(40,16.5,78), Vec(), Vec(), Vec(), Vec(), Refl_t::DIFF)
};




//////////////////////////////////////////////////////////////////////////////
//
//  Class: Plane
//
//  Represents a infinit Plane in 3D space.
//
//////////////////////////////////////////////////////////////////////////////
class Plane : public Object{
private:
	Vec    point;   //a point on the plane;
	Vec    normal;  //normal of plane;

	Vec    color;
	Vec    emission;
	Refl_t reflection;

public:
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



	//constructor
	//Plane()  delete;
	Plane(){};

	Plane( Vec _p, Vec _n, Vec _c, Vec _e,  Refl_t  _reflect){
		this->point      = _p;
		this->normal     = _n;
		this->normal.normalize(); //we want the plane's normal be normalized.

		this->color      = Vec(0.9, 0.0, 0.0);
		//this->color      = _c;
		this->emission   = _e;
		this->reflection = _reflect;
	}

	~Plane(){};


	void get_HitInfo(const Ray &ray, const double t, HitInfo *hit) const{

		hit->object     =  const_cast<Plane *>(this);
		hit->distance   = t;

		hit->position   = ray.o + ray.d*t;
		hit->normal     = this->normal;
		hit->color      = this->color;
		hit->emission   = this->emission;
		hit->reflection = this->reflection;
	}


	double intersect(const Ray &ray) const {
		// n dot (p -p0) = 0;
		// n dot (p' - point) = 0;
		// n dot ( (o + t*d) - point) = 0;
		// n dot (o-p) + t* n dot d = 0;
		// t = (n dot (o -p)) / (n dot d);


		double delta = ray.d.dot(this->normal);
		//double eps = 1e-4;
		double eps = 1;
		if(delta<eps){
			//the ray is parallel to the plane
			return false;
		}
		else{

			Vec op = (ray.o) - (this->point);
			double upper = (this->normal).dot(op);

			double t = upper/delta;

			if(t >= ray.TMax){
				//too far to hit, just return miss
				return 0;
			}
			else{
				return t;
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



Plane planes[]={
  	//Plane( Vec _p,             Vec _n,               Vec _c,                Vec _e,  Refl_t  _reflect);
  	Plane( Vec(50.0, 40.8, 170), Vec(0.0, 0.0, 1.0),   Vec(0.75, 0.75, 0.75), Vec(),  DIFF), //Back
//  
//  	//Ball(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
//  	//Ball(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
//  	//Ball(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back
//  	//Ball(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(.25,.75,.75),DIFF),//Frnt
//  	//Ball(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm
//  	//Ball(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
};




class Ball: public Object{

public:
	double radious;          // radious
	Vec    position;         // position
	Vec    center;           // position
	Vec    emission;         // emission, object self light
	Vec    color;            // color 
	Refl_t reflection;       // reflection type (DIFFuse, SPECular, REFRactive) 


public:

	Ball() = delete;
	Ball(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_): 
		radious(rad_), position(p_), emission(e_), color(c_), reflection(refl_) {};

	~Ball(){}; 

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

		hit->object            = const_cast<Ball *>(this);

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
	
	
//	double intersect(const Ray &r) const { 
//		// returns distance, 0 if nohit 
//		// (p'- p)^2 = R^2;
//		// ((o+t*d) - p)^2 = R^2;
//		// ( t*d + (o-p))^2 - R^2 = 0;
//		// d.d*t^2 + 2*d*(o-p)*t + (o-p).(o-p) - R^2 = 0;
//
//		
//		// Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 
//
//		Vec op    = (this->position) - r.o;       //// direction from light_origin to sphere_center
//		double b  = (op.dot(r.d));
//		double det= b*b-op.dot(op) + radious*radious ; 
//
//		double a1 = r.d.dot(r.d);   //equal 1
//		double b1 = 2*(r.d.dot(r.o-position));
//		double c1 = ((r.o-position).dot(r.o-position)) - radious*radious;
//		double delta = b1*b1 - 4*a1*c1;
//
//
//		if (det<0){
//			assert(delta<0);
//			return 0; 
//		}
//		else{
//			double eps=1e-4;
//			double r0=sqrt(det); 
//			//double t;
//			////return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0); 
//
//			double t0 = b-r0;
//			double t1 = b+r0;
//
//			//to find the positive and minimue one.
//			//obviously t0<t1;
//
//			if(t0>eps){
//				return t0;
//			}
//			else if (t1>eps){
//				return t1;
//			}
//			else{
//				return 0;
//			}
//		}
//	} 


	bool  intersect(const Line &l, Vec3f &enter, Vec3f &exit) const;
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
	
	//Ball(5.0, Vec(40,16.5,78),       Vec(),Vec(1,1,1)*.999, DIFF),//a small ball
//#endif
}; 

class Ellipsoid: public Object{
	//
	//((x*x)/(a*a))+((y*y)/(b*b))+((z*z)/(c*c))=1
	Vec p;
	Vec r;   //(a,b,c)
	
	double intersect(const Ray &ray) const 
	{
			return 0;
	}

	void get_HitInfo(const Ray &ray, const double t, HitInfo *hit) const 
	{

	}
};

//Bezier surface
//Bezier Curve and Surface
class BezierSurface: public Object{
public:
	Vec p;
	Vec r;   //(a,b,c)
	
	double intersect(const Ray &ray) const 
	{
			return 0;
	}

	void get_HitInfo(const Ray &ray, const double t, HitInfo *hit) const 
	{

	}
};

class Curve: public Object{
public:
	Vec p;

	double intersect(const Ray &ray) const 
	{
			return 0;
	}

	void get_HitInfo(const Ray &ray, const double t, HitInfo *hit) const 
	{

	}
};

class Cube: public Object{
private:
	Vec center;       // position
	Vec size;         // size on the length/witdth/height direction

	Vec coords[8];    // Corner coordinates


	Vec emission;     // emission, object self light
	Vec color;        // color 
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive) 

public:
	Cube(){}
	~Cube(){}

	void set_position(){};


	bool intersect(const Ray &ray,  HitInfo* hit) const {
		return false;
	}

	double intersect(const Ray &ray) const 
	{
			return 0;
	}

	void get_HitInfo(const Ray &ray, const double t, HitInfo *hit) const 
	{

	}

};






/////////////////////////////////////////////////////////////


//
//  class Attributes
//  {
//  public:
//  	Attributes() {};
//  	~Attributes() {};
//  
//  
//  	std::string name;
//  	uint32_t start_idx  = 0;      //attri start line index.
//  	uint32_t dimension  = 0;      //vec1? vec2? vec3? vec4?
//  	bool    is_position = false; //whether the attributes is position.
//  	std::vector<glm::vec4> attr;
//  
//  private:
//  
//  };
//




//  class Mesh_config
//  {
//  public:
//  
//  	Mesh_config() {};
//  	~Mesh_config() {};
//  	bool to_draw_mesh_vertex;
//  	bool to_draw_mesh_boundary_vertex;
//  	bool to_draw_mesh_boundary_vertex_connect;
//  
//  	bool to_draw_mesh_line;
//  	bool to_draw_mesh_surface;
//  
//  	bool to_draw_boundingbox_vertex;
//  	bool to_draw_boundingbox_line;
//  	bool to_draw_boundingbox_surface;
//  
//  
//  	//glm::uvec3 mesh_vertex_color            = glm::uvec3(0.0f, 0.0f, 0.0f);// black vertex.
//  	//glm::uvec3 mesh_boundary_vertex         = glm::uvec3(1.0f, 0.0f, 0.0f);//red boundary vertex 
//  	//glm::uvec3 mesh_boundary_vertex_connect = glm::uvec3(43, 243, 66);//green connect
//  
//  	//glm::uvec3 mesh_line                    = glm::uvec3(0.0f, 0.0f, 0.0f); //black line.
//  	//glm::uvec3 mesh_surface                 = glm::uvec3(1.0f, 0.0f, 0.0f);//red surface.
//  
//  	//glm::uvec3 boundingbox_vertex           = glm::uvec3(241, 104, 87);//red boundingbox vertex 
//  	//glm::uvec3 boundingbox_line             = glm::uvec3(0.0f, 1.0f, 0.0f);//green boundingbox line.
//  	//glm::uvec3 boundingbox_surface          = glm::uvec3(1.0f, 0.0f, 0.0f);//red boundingbox surface.
//  
//  	//bool is_need_to_draw[8] = { 0 };//the 8 configure.
//  	//GLuint program[8] = { 0 };//the 8 program.
//  private:
//  
//  };


enum Primitive_Type{
	UNKNOWN      = 0,
	POINT        = 1,
	LINESEGMENT  = 2,
	TRIANGLE     = 3
};

class Mesh{
public:
	//std::string name;

	//Mesh_config config;

	std::string file_name;

	//Generally: we DO NOT consider a mesh is a rigid body.
	//
	bool is_rigid   = false;
	bool is_visible = false;

	//private:
	bool to_get_local_vertex = true;
	uint32_t  local_max_x_vid = -1;
	uint32_t  local_min_x_vid = -1;
	uint32_t  local_max_y_vid = -1;
	uint32_t  local_min_y_vid = -1;
	uint32_t  local_max_z_vid = -1;
	uint32_t  local_min_z_vid = -1;



	uint32_t vertex_count = 0;
	uint32_t mesh_id      = -1;  //openGL EID.

	std::vector<Vec *> vertices;

	Vec * colors   = nullptr;
	Vec * uv       = nullptr;
	
	
	Primitive_Type  primitive_type =  Primitive_Type::UNKNOWN;


	uint32_t primitive_count = 0;
	std::vector<Object *>  primitives;


	void Read_Mesh_Attri();
	void Gen_Object();


	//the buffer only save mesh position.
	//GLfloat *mesh_position_buffer = nullptr;


	static const uint32_t MAX_NUM_OF_ATTRIBUTTES = 16;    // OpenGL spec
	uint32_t num_of_attributes              = 0;                       // how many attributes a mesh have.
	uint32_t num_of_attributes_components   = 0;                       // how many attributes component a mesh have.

	std::vector<std::string> mesh_candidate_in_position_name = { "in_pos","in_Pos", "in_position","in_ATTRIBUTE0", "gl_Position" };

	//std::vector<Attributes> mesh_attributes;
	//Attributes mesh_attributes[MAX_NUM_OF_ATTRIBUTTES];
	//float max_value[3] = { 0.0f, 0.0f, 0.0f };
	//float min_value[3] = { 0.0f, 0.0f, 0.0f };
	//the init model to place the mesh in screen centry.
	//glm::mat4 init_model = glm::mat4(1.0f);



	//glm::vec3  init_trans_vec  = glm::vec3(0.0f); //init shift vec
	//glm::vec3  init_scale_vec  = glm::vec3(1.0f); //init scale vec
	//glm::uvec3 init_rotate_vec = glm::uvec3(0);   //init roate vec

	//the 6 bounding vertex on mesh. which own (max.x, min.x),(max.y, min.y), (max.z, min.z)
	//std::vector<glm::vec4> mesh_boundary_vertex;

	std::vector<uint32_t>  mesh_boundary_vertex_VID;//VID which own max_x, min_x, max_y, min_y, max_z, min_z


	//std::vector<glm::vec4> mesh_boundary_vertex_local;

	//the bounding 8 box vertex
	//std::vector<glm::vec4> mesh_boundingbox_vertex;
	//std::vector<glm::vec4> mesh_boundingbox_vertex_local;

	//PROGRAMS
	//GLuint mesh_vertex_Program                  = 0;
	//GLuint mesh_line_Program                    = 0;
	//GLuint mesh_surface_Program                 = 0;
	//GLuint mesh_boundary_vertex_Program         = 0;
	//GLuint mesh_boundary_vertex_connect_Program = 0;

	//GLuint mesh_boundingBox_vertex_Program      = 0;
	//GLuint mesh_boundingBox_line_Program        = 0;
	//GLuint mesh_boundingBox_surface_Program     = 0;







	//GLuint mesh_surface_VA0                     = 0;
	//GLuint mesh_line_VA0                        = 0;
	//GLuint boundary_vertex_VAO                  = 0;//the 6 boundary vertex on the mesh.

	//GLuint boundingbox_vertex_VA0               = 0; //the 8 boundingbox vertex on the mesh.
	//GLuint boundingbox_surface_VA0              = 0;
	//GLuint boundingbox_line_VA0                 = 0;



	//GLuint Gen_Mesh_Surface_VAO();
	//GLuint Gen_Mesh_Line_VAO();
	//GLuint Gen_Boundary_Vertex_VAO();
	//void Gen_Init_Model_Matrix(bool need_init_model_matrix);
	//
	//
	//void Gen_BoundingBox_Vertex();
	//GLuint Gen_BoundingBox_Vertex_VAO();
	//GLuint Gen_BoundingBox_Surface_VAO();
	//GLuint Gen_BoundingBox_Line_VAO();


	/////////////////////////
//
//#ifdef show_local
//	GLuint boundary_vertex_VAO_local = 0;//the 6 boundary vertex on the mesh.
//
//	GLuint boundingbox_vertex_VAO_local = 0; //the 8 boundingbox vertex on the mesh.
//	GLuint boundingbox_surface_VAO_local = 0;
//	GLuint boundingbox_line_VAO_local = 0;
//
//
//	GLuint Gen_Boundary_Vertex_VAO_local();
//	void Gen_BoundingBox_Vertex_local();
//	GLuint Gen_BoundingBox_Vertex_VAO_local();
//	GLuint Gen_BoundingBox_Surface_VAO_local();
//	GLuint Gen_BoundingBox_Line_VAO_local();
//#endif
//	///////////////////////
//


public:

	Mesh(){};
	~Mesh(){

	};
	Mesh(std::string mesh_input_file_name);
	
	void get_Primitives(std::vector<Object *> &p){
		p=primitives;
	}
};



void Mesh::Read_Mesh_Attri()
{
	auto ClearAllSpace = [&](std::string str) {
		int index = 0;
		if (!str.empty())
		{
			while ((index = str.find(' ', index)) != string::npos)
			{
				str.erase(index, 1);
			}
		}
		return str;
	};


	assert(!file_name.empty() && "error: invalid input mesh file name!");

	uint32_t num_of_token = 0; //each line conatained token number.
	uint32_t line_idx = 0;
	std::string t_line;

	ifstream ifs;
	ifs.open(file_name.data(), std::ifstream::in);
	assert(ifs.is_open() == true && "error: fail to open file.");


	//first line to get split the those attributes.
	if (getline(ifs, t_line))
	{
		line_idx++;
		std::istringstream ss(t_line);


		std::vector<std::string> token_set;
		std::string t_token;
		while (std::getline(ss, t_token, ','))
		{
			token_set.push_back(t_token);
		}


		num_of_token = token_set.size();


		//file header format
		//VTX, IDX, gl_Position.x, gl_Position.y, gl_Position.z, gl_Position.w

		assert( ClearAllSpace(token_set[0]) == "VTX" && "error: fail to parse file header!");
		assert( ClearAllSpace(token_set[1]) == "IDX" && "error: fail to parse file header!");

		assert( ClearAllSpace(token_set[2]) == "gl_Position.x" && "error: fail to parse file header!");
		assert( ClearAllSpace(token_set[3]) == "gl_Position.y" && "error: fail to parse file header!");
		assert( ClearAllSpace(token_set[4]) == "gl_Position.z" && "error: fail to parse file header!");
		assert( ClearAllSpace(token_set[5]) == "gl_Position.w" && "error: fail to parse file header!");
	}
	else
	{
		assert(0 && "fail to get first line.");
	}



	while (getline(ifs, t_line))
	{
		line_idx++;

		std::istringstream ss(t_line);
		std::string token;

		std::vector<string> token_set;
		while (std::getline(ss, token, ',')) {
			token_set.push_back(token);
		}
		assert(token_set.size() == num_of_token);


		//x, y, z, w of a vertices.
		double _x = atof(token_set[2].c_str());
		double _y = atof(token_set[3].c_str());
		double _z = atof(token_set[4].c_str());
		double _w = atof(token_set[5].c_str());

		Vec3f * _v = new Vec3f((_x/_w), (_y/_w), (_z/_w));

		vertices.push_back(_v);
	}

	ifs.close();
	cout << "Read Mesh Attributes success." << endl;

//	std::vector<Vec3f *>::const_iterator ite = vertices.begin();
//
//	if((vertices.size()%3) ==0){
//		//to build triangles
//
//		for( ;ite!= vertices.end(); ite=(ite+3)){
//
//			Vec3f* v0 = (*(ite+0));
//			Vec3f* v1 = (*(ite+1));
//			Vec3f* v2 = (*(ite+2));
//
//			//Triangle * _t = new Triangle((t_element.x/t_element.w), (t_element.y/t_element.w), (t_element.z/t_element.w));
//			//mesh.push_back(_t); 
//		}
//
//	}
//	else{
//		//to build point.
//
//		for( ;ite!= vertices.end(); ite=(ite+1)){
//			Vec3f* v0 = (*(ite+0));
//
//			//Ball * _t = new Ball((t_element.x/t_element.w), (t_element.y/t_element.w), (t_element.z/t_element.w));
//			//mesh.push_back(_t); 
//		}
//	}

	//the 6 value on the mesh.
	float max_x = 0.0f, min_x = 0.0f, max_y = 0.0f, min_y = 0.0f, max_z = 0.0f, min_z = 0.0f;
	uint32_t max_x_vid = 0.0f, min_x_vid = 0.0f, max_y_vid = 0.0f, min_y_vid = 0.0f, max_z_vid = 0.0f, min_z_vid = 0.0f;


	bool is_a_input_mesh  = false;
	bool is_a_output_mesh = false;
	//std::vector<Attributes>::const_iterator attr_ite = mesh_attributes.begin();



//	bool is_find_mesh_pos = false;
//	for (size_t i = 0; i < num_of_attributes; i++)
//	{
//
//		std::string attr_name = mesh_attributes[i].name;
//		if (is_in_position(attr_name) == true) {
//
//			is_find_mesh_pos = true;
//
//			Attributes * attr_ite = &(mesh_attributes[i]);
//			assert((attr_ite->dimension) >= 3 && "error: for any mesh, input position must be greater than 3.");
//
//
//
//			attr_ite->is_position = true;
//
//			// init the 6 mesh boundary vertex
//			max_x = attr_ite->attr[0].x;
//			min_x = attr_ite->attr[0].x;
//
//			max_y = attr_ite->attr[0].y;
//			min_y = attr_ite->attr[0].y;
//
//			max_z = attr_ite->attr[0].z;
//			min_z = attr_ite->attr[0].z;
//
//			for (size_t i = 0; i < 6; i++) {
//				mesh_boundary_vertex.push_back(attr_ite->attr[0]);
//				mesh_boundary_vertex_VID.push_back(-1);
//			}
//
//
//			//to create VAO;
//			//to create VBO;
//			uint32_t t_vertex_count = attr_ite->attr.size();
//
//			//mesh_position_buffer = new GLfloat[t_vertex_count * 4];
//			//memset(mesh_position_buffer, 0, (sizeof(GLfloat) * t_vertex_count * 4));
//
//			std::vector<glm::vec4>::const_iterator ele_ite = attr_ite->attr.begin();
//			for (; ele_ite != attr_ite->attr.end(); ele_ite++) {
//				glm::vec4 t_element = *ele_ite;
//				uint32_t t_vid = ele_ite - attr_ite->attr.begin();
//
//				//mesh_position_buffer[(t_vid * 4) + 0] = t_element.x;
//				//mesh_position_buffer[(t_vid * 4) + 1] = t_element.y;
//				//mesh_position_buffer[(t_vid * 4) + 2] = t_element.z;
//				//mesh_position_buffer[(t_vid * 4) + 3] = t_element.w;
//
//				
//				//OK we find a triangle, save it to the set.
//				Triangle * _t = new Triangle((t_element.x/t_element.w), (t_element.y/t_element.w), (t_element.z/t_element.w));
//				mesh.push_back(_t); 
//
//				//update the boundrary vertex infor
//				if (t_element.x > max_x) {
//					mesh_boundary_vertex[0]     = t_element;
//					mesh_boundary_vertex_VID[0] = t_vid;
//					max_x = t_element.x;
//					max_x_vid = t_vid;
//				}
//
//				if (t_element.x < min_x) {
//					mesh_boundary_vertex[1]     = t_element;
//					mesh_boundary_vertex_VID[1] = t_vid;
//					min_x = t_element.x;
//					min_x_vid = t_vid;
//				}
//
//				if (t_element.y > max_y) {
//					mesh_boundary_vertex[2]     = t_element;
//					mesh_boundary_vertex_VID[2] = t_vid;
//					max_y = t_element.y;
//					max_y_vid = t_vid;
//				}
//
//				if (t_element.y < min_y) {
//					mesh_boundary_vertex[3]     = t_element;
//					mesh_boundary_vertex_VID[3] = t_vid;
//					min_y = t_element.y;
//					min_y_vid = t_vid;
//				}
//
//				if (t_element.z > max_z) {
//					mesh_boundary_vertex[4]     = t_element;
//					mesh_boundary_vertex_VID[4] = t_vid;
//					max_z = t_element.z;
//					max_z_vid = t_vid;
//				}
//
//				if (t_element.z < min_z) {
//					mesh_boundary_vertex[5]     = t_element;
//					mesh_boundary_vertex_VID[5] = t_vid;
//					min_z = t_element.z;
//					min_z_vid = t_vid;
//				}
//
//				//mesh_boundary_vertex[0] = (mesh_boundary_vertex[0].x > t_element.x) ? mesh_boundary_vertex[0] : t_element;
//				//mesh_boundary_vertex[1] = (mesh_boundary_vertex[1].x < t_element.x) ? mesh_boundary_vertex[1] : t_element;
//				//
//				//mesh_boundary_vertex[2] = (mesh_boundary_vertex[2].y > t_element.y) ? mesh_boundary_vertex[2] : t_element;
//				//mesh_boundary_vertex[3] = (mesh_boundary_vertex[3].y < t_element.y) ? mesh_boundary_vertex[3] : t_element;
//				//
//				//mesh_boundary_vertex[4] = (mesh_boundary_vertex[4].z > t_element.z) ? mesh_boundary_vertex[4] : t_element;
//				//mesh_boundary_vertex[5] = (mesh_boundary_vertex[5].z < t_element.z) ? mesh_boundary_vertex[5] : t_element;
//
//
//			}
//			//cout<<"6 boundary value:"<<endl;
//			//cout << max_x << ",\t" << min_x << ",\t" << max_y << ",\t" << min_y << ",\t" << max_z << ",\t" << min_z << endl;
//			//
//
//			/*
//			   static uint32_t mesh_id_1 = 0;
//			   if (mesh_id_1 == 2)
//			   {
//			   for (size_t i = 0; i < 6; i++)
//			   {
//			//init it.
//			mesh_boundary_vertex_local.push_back(attr_ite->attr[0]);
//			}
//			mesh_boundary_vertex_local[0] = attr_ite->attr[local_max_x_vid];
//			mesh_boundary_vertex_local[1] = attr_ite->attr[local_min_x_vid];
//
//			mesh_boundary_vertex_local[2] = attr_ite->attr[local_max_y_vid];
//			mesh_boundary_vertex_local[3] = attr_ite->attr[local_min_y_vid];
//
//			mesh_boundary_vertex_local[4] = attr_ite->attr[local_max_z_vid];
//			mesh_boundary_vertex_local[5] = attr_ite->attr[local_min_z_vid];
//			}
//			mesh_id_1++;
//			 */
//			//
//			//
//			cout << "Mesh boundary vertex:" << endl;
//			cout << "\t max_x_vid: " << max_x_vid << endl;
//			cout << "\t min_x_vid: " << min_x_vid << endl;
//
//			cout << "\t max_y_vid: " << max_y_vid << endl;
//			cout << "\t min_y_vid: " << min_y_vid << endl;
//
//			cout << "\t max_z_vid: " << max_z_vid << endl;
//			cout << "\t min_z_vid: " << min_z_vid << endl;
//
//			cout << "\t max_x_vertex: (" << attr_ite->attr[max_x_vid].x << ", " << attr_ite->attr[max_x_vid].y << ", " << attr_ite->attr[max_x_vid].z << ")" << endl;
//			cout << "\t min_x_vertex: (" << attr_ite->attr[min_x_vid].x << ", " << attr_ite->attr[min_x_vid].y << ", " << attr_ite->attr[min_x_vid].z << ")" << endl;
//
//			cout << "\t max_y_vertex: (" << attr_ite->attr[max_y_vid].x << ", " << attr_ite->attr[max_y_vid].y << ", " << attr_ite->attr[max_y_vid].z << ")" << endl;
//			cout << "\t min_y_vertex: (" << attr_ite->attr[min_y_vid].x << ", " << attr_ite->attr[min_y_vid].y << ", " << attr_ite->attr[min_y_vid].z << ")" << endl;
//
//			cout << "\t max_z_vertex: (" << attr_ite->attr[max_z_vid].x << ", " << attr_ite->attr[max_z_vid].y << ", " << attr_ite->attr[max_z_vid].z << ")" << endl;
//			cout << "\t min_z_vertex: (" << attr_ite->attr[min_z_vid].x << ", " << attr_ite->attr[min_z_vid].y << ", " << attr_ite->attr[min_z_vid].z << ")" << endl;
//
//
//			//max_value[0] = max_x;
//			//max_value[1] = max_y;
//			//max_value[2] = max_z;
//			//
//			//min_value[0] = min_x;
//			//min_value[1] = min_y;
//			//min_value[2] = min_z;
//
//
//			//GLuint attribute_vec_size = 4;
//
//
//			//GLuint  VAO = 0, VBO = 0, EBO = 0;
//			//VERIFY_GL(glGenVertexArrays(1, &VAO));
//			//VERIFY_GL(glGenBuffers(1, &VBO));
//
//
//			//// Bind the Vertex Array Object first, then bind and set vertex buffer(s) and attribute pointer(s).
//			//VERIFY_GL(glBindVertexArray(VAO));
//			//VERIFY_GL(glBindBuffer(GL_ARRAY_BUFFER, VBO));
//			//VERIFY_GL(glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*attribute_vec_size* t_vertex_count, mesh_position_buffer, GL_STATIC_DRAW));
//
//			//GLint pos_location = -1;
//			//VERIFY_GL_RET(pos_location, glGetAttribLocation(mesh_surface_Program, "in_position"));
//			//VERIFY_GL(glEnableVertexAttribArray(pos_location));
//			//VERIFY_GL(glVertexAttribPointer(pos_location, attribute_vec_size, GL_FLOAT, GL_FALSE, 0, (GLvoid*)0));
//			//VERIFY_GL(glDisableVertexAttribArray(pos_location));
//
//			//VERIFY_GL(glBindBuffer(GL_ARRAY_BUFFER, 0)); // Note that this is allowed, the call to glVertexAttribPointer registered VBO as the currently bound vertex buffer object so afterwards we can safely unbind
//
//			///*
//			//   VERIFY_GL(glGenBuffers(1, &EBO));
//			//   VERIFY_GL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO));
//			//   VERIFY_GL(glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW));
//			//   VERIFY_GL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));
//			// */
//
//			//VERIFY_GL(glBindVertexArray(0)); // Unbind VAO (it's always a good thing to unbind any buffer/array to prevent strange bugs), remember: do NOT unbind the EBO, keep it bound to this VAO
//
//
//
//
//
//			//cout << "Gne Mesh Surface VAO success." << endl;
//			//cout << "\t Mesh boundary vertex is:" << endl;
//			//cout << "\t \t max_x vertex is :" << "(" << mesh_boundary_vertex[0].x << ", " << mesh_boundary_vertex[0].y << "," << mesh_boundary_vertex[0].z << ")" << endl;
//			//cout << "\t \t min_x vertex is :" << "(" << mesh_boundary_vertex[1].x << ", " << mesh_boundary_vertex[1].y << "," << mesh_boundary_vertex[1].z << ")" << endl;
//			//cout << endl;
//			//cout << "\t \t max_y vertex is :" << "(" << mesh_boundary_vertex[2].x << ", " << mesh_boundary_vertex[2].y << "," << mesh_boundary_vertex[2].z << ")" << endl;
//			//cout << "\t \t min_y vertex is :" << "(" << mesh_boundary_vertex[3].x << ", " << mesh_boundary_vertex[3].y << "," << mesh_boundary_vertex[3].z << ")" << endl;
//			//cout << endl;
//			//cout << "\t \t max_z vertex is :" << "(" << mesh_boundary_vertex[4].x << ", " << mesh_boundary_vertex[4].y << "," << mesh_boundary_vertex[4].z << ")" << endl;
//			//cout << "\t \t min_z vertex is :" << "(" << mesh_boundary_vertex[5].x << ", " << mesh_boundary_vertex[5].y << "," << mesh_boundary_vertex[5].z << ")" << endl;
//			//cout << endl;
//
//			//return VAO;
//		}
//	}
//	assert(is_find_mesh_pos  != false && "error: fail find mesh position.");

}



void Mesh::Gen_Object(){

	vertex_count = vertices.size();


	if(vertex_count%3==0){
		primitive_type = Primitive_Type::TRIANGLE;
	}
	else if (vertex_count%2==0){ 
		primitive_type = Primitive_Type::LINESEGMENT;
	}
	else{
		assert(0 && "error: not suport.");
		primitive_type = Primitive_Type::POINT;
	}


	switch(primitive_type){

		case Primitive_Type::POINT:
			{
				primitive_count = vertex_count;

				for(uint32_t i=0; i < primitive_count; i++){
					Point * p = new Point();
					//vertices[i], 
					primitives.push_back(p);
				}

			}
			break;

		case Primitive_Type::LINESEGMENT:
			{
				assert(vertex_count % 2 ==0 && "error: for triangles, the vertex count must eqaul 2*n.");
				primitive_count = (vertex_count/2);

				for(uint32_t i=0; i <primitive_count; i++){
					//LineSegment(const Vec3f &_start, const Vec3f &_end, const Vec3f &_color, double _thickness);

					//default is black linesegment.
					LineSegment * l = new LineSegment( (*vertices[(i*2)+0]),  (*vertices[(i*2)+1]),  Vec3f(),  0.1);
					primitives.push_back(l);
				}
			}
			break;

		case Primitive_Type::TRIANGLE:
			{
				assert(vertex_count % 3 ==0 && "error: for triangles, the vertex count must eqaul 3*n.");

				primitive_count = (vertex_count/3);

				for(uint32_t i=0; i < primitive_count; i++){
					//vertices[i*3+0]; 
					//vertices[i*3+1];
					//vertices[i*3+2];

					Vec _c[3]={Vec(1,1,1)*.999, Vec(1,1,1)*.999, Vec(1,1,1)*.999};
					Vec _e;

					Refl_t _r = Refl_t::DIFF;

					Triangle * t = new Triangle( (*vertices[(i*3)+0]),  (*vertices[(i*3)+1]), (*vertices[(i*3)+2]),  _c[0], _c[1], _c[2], _e, _r);
					primitives.push_back(t);
				}
			}
			break;


		default:
			assert(0 && "error: unknow primtive type.");
	}

	//remove all sapce in string.
//	auto ClearAllSpace = [&](std::string str) {
//		int index = 0;
//		if (!str.empty())
//		{
//			while ((index = str.find(' ', index)) != string::npos)
//			{
//				str.erase(index, 1);
//			}
//		}
//		return str;
//	};



//	auto is_in_position = [&](std::string name) {
//
//		bool res = false;
//		std::string clean_name = ClearAllSpace(name);
//
//		uint32_t candiate_number = mesh_candidate_in_position_name.size();
//		std::vector<std::string>::const_iterator ite = mesh_candidate_in_position_name.begin();
//		for (; ite != mesh_candidate_in_position_name.end(); ite++)
//		{
//
//			if (clean_name == *ite)
//			{
//				res = true;
//				return res;
//			}
//		}
//		return res;
//	};
//
//
//	//the 6 value on the mesh.
//	float max_x = 0.0f, min_x = 0.0f, max_y = 0.0f, min_y = 0.0f, max_z = 0.0f, min_z = 0.0f;
//	uint32_t max_x_vid = 0.0f, min_x_vid = 0.0f, max_y_vid = 0.0f, min_y_vid = 0.0f, max_z_vid = 0.0f, min_z_vid = 0.0f;
//
//
//	bool is_a_input_mesh = false;
//	bool is_a_output_mesh = false;
//	//std::vector<Attributes>::const_iterator attr_ite = mesh_attributes.begin();



//	bool is_find_mesh_pos = false;
//	for (size_t i = 0; i < num_of_attributes; i++)
//	{
//
//		std::string attr_name = mesh_attributes[i].name;
//		if (is_in_position(attr_name) == true) {
//
//			is_find_mesh_pos = true;
//
//			Attributes * attr_ite = &(mesh_attributes[i]);
//			assert((attr_ite->dimension) >= 3 && "error: for any mesh, input position must be greater than 3.");
//
//
//
//			attr_ite->is_position = true;
//
//			// init the 6 mesh boundary vertex
//			max_x = attr_ite->attr[0].x;
//			min_x = attr_ite->attr[0].x;
//
//			max_y = attr_ite->attr[0].y;
//			min_y = attr_ite->attr[0].y;
//
//			max_z = attr_ite->attr[0].z;
//			min_z = attr_ite->attr[0].z;
//
//			for (size_t i = 0; i < 6; i++) {
//				mesh_boundary_vertex.push_back(attr_ite->attr[0]);
//				mesh_boundary_vertex_VID.push_back(-1);
//			}
//
//
//			//to create VAO;
//			//to create VBO;
//			uint32_t t_vertex_count = attr_ite->attr.size();
