

ray{
	vec3 o;
	vec3 d;
}


ray = o+ t*d;

sphere:{
 vec3 c;

 double r;

 (p-c)^2 - r^2= 0;
}

f(o+td)=o;

//三角形求教

// triangle & ray intersaction.

//geometry inside/outside

//ray intersaction will surface, point in 

//重心坐标，MT 方法，克莱姆法则
//在三角形内，

//accelerate interacater way.

AABB

AIX Algin BoundingBox (AABB)

//线段求交集

tenter =max{tmin}, texit = min(tmax)

ray(ray){

	if(ray hit on any object){
	}

}



intersection(ray r, BVH node){
	if(r mis on node.bbaa){
		return;
	}
	else{
		if(node is leaf){
			for(all objects in node){
				find the closet object.
			}
		}
		else{

			hit1 = intersection(r, node.child1);
			hit1 = intersection(r, node.child2);

			return hit closet of hit1, hit2
		}
	}
}
