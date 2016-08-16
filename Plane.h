#ifndef __PLANE_H__
#define __PLANE_H__

#include "CKAll.h"
#include "MyUtility.h"

//http://www.cnblogs.com/graphics/archive/2009/10/17/1585281.html
//http://www.tuicool.com/articles/rUrMJvi

//a point and a norm determines a plane

enum Direction{
	Outside,
	Inside
};

class Plane
{
public:
	Plane(){
	}

	Plane(VxVector _norm,VxVector _point):norm(_norm),point(_point){
	}

	bool LineIntersectPlane(VxVector pt,VxVector line,VxVector& ret){
		float dotP = VxVectorInnerProduct(line,norm);
		if(dotP == 0.0f){
			return false;
		}
		else{
			float t = VxVectorInnerProduct(point - pt,norm)/dotP;
			ret = pt + (t * line);
			return true;
		}
	}

	bool RayIntersectTest(VxVector pt,VxVector line){
		float dotP = VxVectorInnerProduct(line,norm);
		if(dotP == 0.0f ){
			return false;
		}
		else{
			float t = VxVectorInnerProduct(point - pt,norm)/dotP;
			if(t>=0)
			    return true;
			else
				return false;
		}
	}

	bool RayIntersectTest(VxVector pt,VxVector line,VxVector& ret){
		float dotP = VxVectorInnerProduct(line,norm);
		if(dotP == 0.0f){
			/*
			float d = - VxVectorInnerProduct(norm , point);
			float dis = VxVectorInnerProduct(norm , pt) + d;
			float direction = 
			if(dis < 0 - 0.001){
				throw "dis < 0 - 0.001";
			    return false;
			}else{
				ret = pt - dis * norm;
				return true;
			}
			*/
			return false;
		}
		else{
			float t = VxVectorInnerProduct(point - pt,norm)/dotP;
			if(t> 0){
				ret = pt + (t * line);
				return true;
			}
			else{
				ret = pt + (t * line);
				return false;
			}
		}
	}

	bool RayIntersectTest(VxVector pt,VxVector line,bool& onThePlane,bool& isParallel,VxVector& ret){
		float dotP = VxVectorInnerProduct(line,norm);
		if(fabs(dotP) < 0.01){
			//throw std::string("RayIntersectTest  ") + std::string(__FILE__);
			//onThePlane = true;
			if(fabs(VxVectorInnerProduct(point - pt,norm)) < 0.01){
			    isParallel = true;
				ret = pt;
				return true;
			}
			return false;
		}
		else{
			float t = VxVectorInnerProduct(point - pt,norm)/dotP;
			if(fabs(t) < 0.01){
				onThePlane = true;
				ret = pt + (t * line);
			    return true;
			}else if(t> 0){
				ret = pt + (t * line);
				return true;
			}
			else{
				ret = pt + (t * line);
				return false;
			}
		}
	}


	Direction DirectionOfLineToPlane(VxVector line){
		float Cos = GetCosAngleOfVectors(line,norm);
		if(Cos > 0)
			return Outside;
		else
			return Inside;
	}

	VxVector norm;
	VxVector point;
};

#endif