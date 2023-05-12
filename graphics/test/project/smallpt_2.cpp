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
//			GLuint attribute_vec_size = 4;
//
//
//			GLuint  VAO = 0, VBO = 0, EBO = 0;
//			VERIFY_GL(glGenVertexArrays(1, &VAO));
//			VERIFY_GL(glGenBuffers(1, &VBO));
//
//
//			// Bind the Vertex Array Object first, then bind and set vertex buffer(s) and attribute pointer(s).
//			VERIFY_GL(glBindVertexArray(VAO));
//			VERIFY_GL(glBindBuffer(GL_ARRAY_BUFFER, VBO));
//			VERIFY_GL(glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*attribute_vec_size* t_vertex_count, mesh_position_buffer, GL_STATIC_DRAW));
//
//			GLint pos_location = -1;
//			VERIFY_GL_RET(pos_location, glGetAttribLocation(mesh_surface_Program, "in_position"));
//			VERIFY_GL(glEnableVertexAttribArray(pos_location));
//			VERIFY_GL(glVertexAttribPointer(pos_location, attribute_vec_size, GL_FLOAT, GL_FALSE, 0, (GLvoid*)0));
//			VERIFY_GL(glDisableVertexAttribArray(pos_location));
//
//			VERIFY_GL(glBindBuffer(GL_ARRAY_BUFFER, 0)); // Note that this is allowed, the call to glVertexAttribPointer registered VBO as the currently bound vertex buffer object so afterwards we can safely unbind
//
//			/*
//			   VERIFY_GL(glGenBuffers(1, &EBO));
//			   VERIFY_GL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO));
//			   VERIFY_GL(glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW));
//			   VERIFY_GL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));
//			 */
//
//			VERIFY_GL(glBindVertexArray(0)); // Unbind VAO (it's always a good thing to unbind any buffer/array to prevent strange bugs), remember: do NOT unbind the EBO, keep it bound to this VAO
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
//			return VAO;
//		}
//	}
//	assert(is_find_mesh_pos  != false && "error: fail find mesh position.");
//	assert(mesh_vertex_count != 0     && "error: fail to reconginize the mesh.");
//  return 0;

}



Mesh::Mesh(std::string mesh_file_name){
	cout << "Creating Mesh:" << mesh_file_name << std::endl;
	file_name = mesh_file_name;

	Read_Mesh_Attri();
	Gen_Object();
}



void Load_Mesh(std::string path, std::string file_suffix, std::vector<Mesh *> &meshes){

	std::vector<std::string> files;
	struct dirent *entry;
	DIR *dir = opendir(path.c_str());

	if (dir == NULL)
	{
		//return files;
		assert(0 && "fail to open dir.");
	}

	while ((entry = readdir(dir)) != NULL)
	{
		files.push_back(entry->d_name);

		std::string file_name = entry->d_name;

		if(file_name.find(file_suffix)!=std::string::npos){
			std::string mesh_file_name = path+"/"+file_name;
			Mesh *m = new Mesh(mesh_file_name);
			meshes.push_back(m);
		}

	}
	closedir(dir);

	//for(std::vector<std::string>::iterator ite=files.begin(); ite!=files.end();ite++){
	//}

}







/////////////////////////////////////////////////////////////
class BVH_Node{
private:
	// AABB: aix-align-boundingbox
	// 0~5: min_x, max_x, min_y, max_y, min_z, max_z
	double aabb[6]={0};

	//two child
	BVH_Node *child0 = nullptr;
	BVH_Node *child1 = nullptr;



	//objects inside the node
	uint32_t obj_cnt = 0;
	Object * obj_set[MAX_OBJECT_CNT_IN_A_BVH_LEAF];


private:
	//__attribute((always_inline)) 
	inline bool is_intersect_aabb(const Ray& ray){

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
	
	
	






public:
	BVH_Node(){
		child0 = nullptr;
		child1 = nullptr;
		memset(aabb,0,sizeof(double)*6);
		
		obj_cnt = 0;
		memset(obj_set, 0, sizeof(Object *)*MAX_OBJECT_CNT_IN_A_BVH_LEAF);
	};

	~BVH_Node(){
		if(!child0){
			delete child0; 
			child0 =nullptr;
		}

		if(!child1){
			delete child1; 
			child1 =nullptr;
		}
		
		memset(aabb,0,sizeof(double)*6);

		obj_cnt = 0;
		memset(obj_set, 0, sizeof(Object *)*MAX_OBJECT_CNT_IN_A_BVH_LEAF);
	};

	//init aabb to the bvh node
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

	//add an object to a leaf_node.
	void add_object(Object * obj){	

		assert((child0 == nullptr) && "error: can't add object to mid_node!");
		assert((child1 == nullptr) && "error: can't add object to mid_node!");
		
		assert(obj_cnt<=MAX_OBJECT_CNT_IN_A_BVH_LEAF && "error: to much objects in the leaf_node!" );

		obj_set[obj_cnt] = obj;
		obj_cnt++;
	}


	//reset object set to a leaf_node.
	void set_object(const std::vector<Object *> &src_set){

		assert((child0 == nullptr) && "error: can't add object to mid_node!");
		assert((child1 == nullptr) && "error: can't add object to mid_node!");

		uint32_t cnt = src_set.size();

		for(uint32_t i=0; i<cnt;i++){

			assert(obj_cnt<=MAX_OBJECT_CNT_IN_A_BVH_LEAF && "error: to much objects in the leaf_node!" );
			
			Object * obj_ptr = src_set[i];
			
			obj_set[obj_cnt] = obj_ptr;
			obj_cnt++;
		}
	}


	void set_child(BVH_Node *p_child0, BVH_Node *p_child1){
		child0 = p_child0;
		child1 = p_child1;
	}


	//do ray bvh intersection.

	bool intersect(const Ray &ray, double &nearest_t, Object *&nearest_obj){
		
	
		if (!is_intersect_aabb(ray)){
			//miss on the node's AABB, return miss.
			return false;
		}
		else{
			//hit on the node's AABB, continue check whether it hit on any object

			if( (this->child0 != nullptr )){
				// this is a mid_node, for any mid_node it must have two child.

				// continue traversal all child node, to find the nearest hit point.
				double t0, t1;
				Object *obj0, *obj1;

				bool res0 = child0->intersect(ray, t0, obj0);
				bool res1 = child1->intersect(ray, t1, obj1);

				if((res0==false) && (res1==false)){
					//miss on two child node
					return false;
				}
				else if((res0==true) && (res1==false)){
					//hit on child0;
					nearest_t   = t0;
					nearest_obj = obj0;
					return true;
				}
				else if((res0==false) && (res1==true)){
					//hit on child1;
					nearest_t   = t1;
					nearest_obj = obj1;
					return true;
				}
				else if((res0==true) && (res1==true)){
					//two child both hit, to find the nearest one.

					if(t0 < t1){
						//hit on child0;
						nearest_t   = t0;
						nearest_obj = obj0;
					}
					else{
						//hit on hit1;
						nearest_t   = t1;
						nearest_obj = obj1;
					}
					return true;
				}

			}
			else{

				//this is a BVH leaf_node 
				//to test all objects in the leaf node, and find the nearest one.
				nearest_obj = nullptr;
				nearest_t   = MAX_RAY_DISTANCE;
				double inf  = MAX_RAY_DISTANCE;

				for(uint32_t i=0;i<obj_cnt;i++){

					double t = obj_set[i]->intersect(ray);
					
					// hit on a object, we need find the nearest one.
					if( (t!=0) && (t< nearest_t) ){
						//update the nearest t and object
						nearest_t   = t;
						nearest_obj = obj_set[i];
					}
				}

				return (nearest_t < inf);



				//	ray.dump();

				//	if(hit_obj.object == nullptr){
				//		//the ray DO NOT hit on any object in the BVH node, return miss

				//		if(b_debug){
				//			log_ofs<<"\t Miss. \n:"<<endl;
				//		}
				//		return false;
				//	}
				//	else{

				//		if(b_debug){

				//			if(dynamic_cast<Ball*>(hit_obj.object)){
				//				log_ofs<<"\t Hit on a Ball:"<<endl;
				//			}

				//			if(dynamic_cast<Plane*>(hit_obj.object)){
				//				log_ofs<<"\t Hit on a Plane:"<<endl;
				//			}

				//			if(dynamic_cast<Triangle*>(hit_obj.object)){
				//				log_ofs<<"\t Hit on a Triangle:"<<endl;
				//			}

				//			log_ofs.flush();
				//			hit_obj.dump();
				//		}

				//		*hit = hit_obj;
				//		return true;

				//	}

			}
		}

	}
};

//////////////////////////////////////////////////////////////////////////////////////////
bool compare_on_x(Object * i1, Object * i2){
	return ((i1->max_x()) < (i2->max_x()));
}

bool compare_on_y(Object * i1, Object * i2){
	return ((i1->max_y()) < (i2->max_y()));
}

bool compare_on_z(Object * i1, Object * i2){
	return ((i1->max_z()) < (i2->max_z()));
}

class Scene{
private:
	uint32_t idx;
	
	//all objects need contatined in the set.
	std::vector<Object *> objects;

	//all inifit objects need contatined in the set.
	std::vector<Object *> infinit_objects;



	//all meshes saved in the set
	std::set<Mesh *> meshes;


	//BVH root ptr.
	BVH_Node *bvh_root = nullptr;

	//all leaf need save here, 
	//some object need add to leaf directly.
	std::vector<BVH_Node *> leaf_nodes;


private:
	/////////////////////////////////////////////////////

	void find_aabb(const std::vector<Object *> &object_set, double *_t_aabb){

		//to find the global longest axis.
		for(size_t i=0; i<object_set.size(); i++){

			Object * obj = object_set[i];

			double t[6] = {0};//boundray of current object.

			t[0] = obj->min_x();
			t[1] = obj->max_x();

			t[2] = obj->min_y();
			t[3] = obj->max_y();

			t[4] = obj->min_z();
			t[5] = obj->max_z();

			_t_aabb[0] = (_t_aabb[0]<t[0])?_t_aabb[0]:t[0];//min_x
			_t_aabb[1] = (_t_aabb[1]>t[1])?_t_aabb[1]:t[1];//max_x

			_t_aabb[2] = (_t_aabb[2]<t[2])?_t_aabb[2]:t[2];//min_y
			_t_aabb[3] = (_t_aabb[3]>t[3])?_t_aabb[3]:t[3];//max_y

			_t_aabb[4] = (_t_aabb[4]<t[4])?_t_aabb[4]:t[4];//min_z
			_t_aabb[5] = (_t_aabb[5]>t[5])?_t_aabb[5]:t[5];//max_z
		}
	}


	void split_object(
			std::vector<Object *> &src_set, 
			std::vector<Object *> &set0,
			std::vector<Object *> &set1,
			double *_aabb
			){

		size_t  cnt = src_set.size();
		assert( cnt!=0 && "invalid src_set");


		//greater than 10 continue split.
		double x_len = (_aabb[1] - _aabb[0]);
		double y_len = (_aabb[3] - _aabb[2]);
		double z_len = (_aabb[5] - _aabb[4]);


		if( (x_len >= y_len) && (x_len >= z_len) ){
			//sort on X axis
			std::sort(src_set.begin(), src_set.end(), compare_on_x);
		}
		else if( (y_len >= x_len) && (y_len >= z_len) ){
			//sort on Y axis
			std::sort(src_set.begin(), src_set.end(), compare_on_y);
		}
		else {
			//sort on Z axis
			std::sort(src_set.begin(), src_set.end(), compare_on_z);
		}


		//split the object to two subset
		size_t half = cnt/2;
		for(size_t i=0;i<cnt;i++){

			Object *obj = src_set[i];

			if(i <= half){
				set0.push_back(obj);
			}
			else{
				set1.push_back(obj);
			}
		}
	}
	








	//to create bvh tree, return root ptr.
	BVH_Node *create_bvh_tree_helper(std::vector<Object *> &object_set){

		size_t count = object_set.size();
		assert(count!=0 && "error: invalid object count!");

		//create a BVH_Node for those object.
		//NOTE: Here must use the "new" operator, DO NOT use local variable.
		BVH_Node * p_node = new BVH_Node();


		// 6 boundary value of all those objects on three axis
		//AABB aabb = new AABB( min_x, max_x, min_y, max_y, min_z, max_z);
		double t_aabb[6]={
			DBL_MAX, DBL_MIN, 
			DBL_MAX, DBL_MIN, 
			DBL_MAX, DBL_MIN};

		//aabb value saved in the buffer
		find_aabb(object_set, t_aabb);

		//init aabb
		p_node->init_aabb(t_aabb);


		
		if(count> MAX_OBJECT_CNT_IN_A_BVH_LEAF){
			//object more than threshold, continue split objects.


			std::vector<Object *> set0, set1;
			split_object(object_set, set0, set1, t_aabb);

			BVH_Node * p_child0 = create_bvh_tree_helper(set0);
			BVH_Node * p_child1 = create_bvh_tree_helper(set1);

			//set the children ptr to the node.
			p_node->set_child(p_child0, p_child1);
			
			cout<<"\t MID_BVH_NODE: "<< p_node <<" , object count = "<< count<<endl;
		}
		else{
			//less than threshold, this node is leaf
			p_node->set_object(object_set);

			//we need to know all leaf_node
			leaf_nodes.push_back(p_node);

			cout<<"\t BVH_LEAF: "<< p_node <<" , object count = "<< count<<endl;
		}

		return p_node;
	}

	void gen_bvh_tree(){

		//build bvh gen BVH_TREE
		if(bvh_root!=nullptr){

			delete bvh_root;
			bvh_root = nullptr;
		}
		bvh_root = create_bvh_tree_helper(objects);
	}

public:

	void clear(){
		if(bvh_root!=nullptr){
			delete bvh_root;
			bvh_root = nullptr;
		}
		
		leaf_nodes.clear();
		meshes.clear();
		objects.clear();
	}


	Scene(){
		bvh_root = nullptr;
		
		leaf_nodes.clear();
		meshes.clear();
		objects.clear();
	}

	~Scene(){
		clear();
	}


	void add_Objects(const std::vector<Object *> &objs){

		for(std::vector<Object *>::const_iterator ite = objs.begin(); ite!=objs.end(); ite++){
			Object * _p = *ite;
			objects.push_back(_p);
		}
	}

	void add_Mesh(Mesh *m){

		assert(m!=nullptr && "error: invalid mesh ptr.");

		std::set<Mesh *>::iterator ite = meshes.find(m);
		if(ite==meshes.end()){
			meshes.insert(m);
		}
	}


	void delete_Mesh(Mesh *m){
		assert(m!=nullptr && "error: invalid mesh ptr.");

		std::set<Mesh *>::iterator ite = meshes.find(m);
		if(ite!=meshes.end()){
			meshes.erase(m);
		}
	}

	void add_Meshes(std::vector<Mesh *> &_meshes){
		for( std::vector<Mesh *>::iterator ite = _meshes.begin(); ite!=_meshes.end();ite++){
			Mesh *m = (*ite);
			add_Mesh(m);
		}
	}

	void delete_Meshes(std::vector<Mesh *> &_meshes){
		for( std::vector<Mesh *>::iterator ite = _meshes.begin(); ite!=_meshes.end();ite++){
			Mesh *m = (*ite);
			delete_Mesh(m);
		}
	}




	void add_Mesh_Primitives(){
		for(std::set<Mesh *>::iterator ite=meshes.begin(); ite!=meshes.end(); ++ite){
			Mesh * m = *ite;

			std::vector<Object *>  p;
			m->get_Primitives(p);

			//add mesh to the Object set.
			add_Objects(p);
		}
	}





	//return ray_color
	Vec radiance_1(const Ray &ray,  unsigned short *Xi){ 

		double t             = DBL_MAX;
		Object * nearest_obj = nullptr;
		if(!(bvh_root->intersect(ray, t, nearest_obj))){
			// the ray do not hit on any object, just return black.
			return Vec();
		}
		else{

			// the ray hit on a object,
			// continue trace the ray OR return the emission of object.

			HitInfo hitinfo;
			nearest_obj->get_HitInfo(ray, t, &hitinfo);

			Vec f  = hitinfo.color;                  // object color, (BRDF modulator) 

			//max refl 
			//max color-component on the oject.
			double p = (f.x>f.y && f.x>f.z) ? f.x : ((f.y>f.z) ? f.y : f.z); 


			if (ray.depth()>(MAX_RAY_FELECT_CNT-1)){ 

				//	//reflect too much times, must return now.
				//	if(ray.depth()>MAX_RAY_FELECT_CNT_MUST_RETURN){

				//		//return obj.e; 
				//		//NOTE: Hit on the object, return object emission.
				//		return hit_obj.emission;
				//	}
				//	else{
				//	}

				//Russian Roulette, RR 
				if(erand48(Xi)>p){ 
					//return obj.e; 
					//NOTE: Hit on the object, return object emission.
					return hitinfo.emission;

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

			switch(hitinfo.reflection){

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

						//Vec w = hit_obj.normal;                                               //z direction
						//Vec u =((fabs(w.x)>(0.1)? Vec(0,1,0): Vec(1,0,0)).cross(w)).norm();   //u is perpendicular to w
						//Vec v = w.cross(u);    //cross.                                       //v is perpendicular to u, w
						////random Sampling Unit Hemisphere 1-r2;
						//Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); // d is a random reflect direction.

						////Vec w=nl;
						Vec w = hitinfo.normal;                                               //z direction
						Vec u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(); 
						Vec v= w%u; 
						//Vec v = w.cross(u);    //cross.                                       //v is perpendicular to u, w
						Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); 




						//get the diffuse_ray
						//NOTE: the random ray's depth need "+1".
						uint32_t new_depth = ray.depth() + 1; 
						Ray diffuse_ray = Ray(hitinfo.position, d, new_depth);

						//step 2: Sampling Unit Hemisphere
						//recursively
						Vec diffuse_color =  this->radiance_1(diffuse_ray, Xi);

						return (hitinfo.emission + f.mult(diffuse_color));



						// Ideal DIFFUSE reflection
						//  double r1=2*M_PI*erand48(Xi);
						//  double r2=erand48(Xi);
						//  double r2s=sqrt(r2);
						//  Vec w=nl;
						//  Vec u=((fabs(w.x)>0.1?Vec(0,1,0):Vec(1,0,0))%w).norm();
						//  Vec v= w%u;
						//  Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
						//  return obj.e + f.mult(radiance(Ray(x,d),depth,Xi));

					}
					break;

				case SPEC:
					{
						// Ideal SPECULAR reflection
						assert(0 && "never go here");




						//step 1. get a reflect ray
						//NOTE: the reflection ray's depth need "+1".
						// Ray reflect_ray = ray.reflect(hit_obj.position, hit_obj.normal);

						Vec reflect_dir;//FIXME: todo
						reflect_dir.norm();

						uint32_t new_depth = ray.depth() + 1;  
						Ray reflect_ray = Ray(hitinfo.position, reflect_dir, new_depth);

						//step 2. continue trace the reflect_ray recursively
						Vec spec_color =  this->radiance_1(reflect_ray, Xi);

						return hitinfo.emission + f.mult(spec_color);
					}
					break;

				case REFR:
					{

						// Ideal dielectric REFRACTION
						assert(0 && "un-supported.");

#if 0
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


	//render a camera's image
	void render(){
		//
		//
		//
	};



	//the main function to render a frame
	void render_frame_1(
		Vec cam_pos /*camera pos*/,
		Vec cam_dir /*camera direction*/,
		Vec *c     /*buffer to save frame*/, 

		const uint32_t w /*image width*/, 
		const uint32_t h /*image hight*/,
		const uint32_t samps /*samples per subpixel*/
		){

		add_Mesh_Primitives();
		gen_bvh_tree();
		
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

							//a new ray, depth is 0.
							Ray ray = Ray(cam.o+d*140, d.norm(), 0);


							//we get the ray's (r,g,b)
							//Vec t = radiance(Ray(cam.o+d*140,d.norm()),0,Xi); 
							Vec t = radiance_1(ray, Xi); 

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
};























int main(int argc, char *argv[]){ 


	uint32_t w=1024; 
	uint32_t h=768;
	uint32_t samps = ((argc==2) ? atoi(argv[1])/4 : 1); // # samples 


	double delta_t = ((double)(1.0))/FPS;

	uint32_t frame_cnt = total_t/delta_t; //how many frame we need to render.

	
	Vec *c=new Vec[w*h]; //buffer to save finnal image. the init value is (0,0,0)
	
	Vec init_pos = spheres[9].p;
	Vec center   = spheres[7].p;

	uint32_t frame_idx = 0;
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
		for(uint32_t idx=0;idx< (sizeof(balls)/sizeof(Ball)); idx++){
			src_set.push_back( &(balls[idx]) );
		}

		for(uint32_t idx=0; idx< sizeof(triangles)/sizeof(Triangle); idx++){
			src_set.push_back( &(triangles[idx]) );
		}

	//	for(uint32_t idx=0; idx< sizeof(planes)/sizeof(Plane); idx++){
	//		src_set.push_back( &(planes[idx]) );
	//	}

		std::string dir="./F3";
		//std::string dir="./test";
		std::string file_suffix = "_out.csv";

		std::vector<Mesh *> meshes;
		Load_Mesh(dir,file_suffix, meshes);


		Scene * t_scene = new Scene();


		t_scene->add_Objects(src_set);

		t_scene->add_Meshes(meshes);
		
		//render
		t_scene->render_frame_1(cam_pos, cam_dir, c,  w, h, samps);

		//render
		//render_frame_0(cam_pos, cam_dir, c, (&bvh), w, h, samps);

		//save
		save_frame(c, w, h, frame_idx);
	}
	
	//free buffer.
	delete []c;

	//generage a gif
	convert_ppm_to_gif();

	return 0;
