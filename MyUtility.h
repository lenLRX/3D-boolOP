#ifndef __MYUTILITY_H__
#define __MYUTILITY_H__

#include "vxvector.h"
#include <string>

static float VxVectorInnerProduct(VxVector v1,VxVector v2){
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

static VxVector VxVectorCrossProduct(VxVector a, VxVector b){
	return VxVector(a.y * b.z - a.z *b.y ,
		            a.z * b.x - a.x * b.z,
					a.x * b.y - a.y *b.x);
}

//Angle in Rad, not Degree!!
static float GetCosAngleOfVectors(VxVector v1, VxVector v2){
	float ProductOfMagnitude = v1.Magnitude()*v2.Magnitude();
	if(ProductOfMagnitude != ProductOfMagnitude)
		return ProductOfMagnitude;//return NaN to caller ,let Them handle it,I really want to throw a Exception!

	return VxVectorInnerProduct(v1,v2)/ProductOfMagnitude;
}


struct VxVectorLess{
	bool operator ()(const VxVector& lhs,const VxVector& rhs)const {
		if(lhs.x < rhs.x){
			return true;
		}else if(lhs.x == rhs.x){
			if(lhs.y < rhs.y){
				return true;
			}else if(lhs.y == rhs.y){
				if(lhs.z < rhs.z)
					return true;
				else
					return false;
			}else
				return false;
		}else
			return false;
	}
};

static bool AlmostEqualZero(float z){
	return fabs(z)<0.001;
}

static bool AlmostEqual(float lhs,float rhs){
	return fabs(lhs-rhs) <0.001;
	/*
	if(lhs == 0.0f){
		if(rhs == 0.0f)
			return true;
		else
			return fabs((lhs-rhs)/rhs) <0.001;
	}else
	    return fabs((lhs-rhs)/lhs) <0.001;
    */
}

static bool SameVertex(VxVector v1,VxVector v2){
	return AlmostEqual(v1.v[0],v2.v[0]) && AlmostEqual(v1.v[1],v2.v[1]) && AlmostEqual(v1.v[2],v2.v[2]);
}

struct FaceTuple
{
	FaceTuple(CKVINDEX f1,CKVINDEX f2,CKVINDEX f3){
		fIndex[0] = f1;
		fIndex[1] = f2;
		fIndex[2] = f3;
	}

	CKVINDEX fIndex[3];
};

static float distance2line(VxVector a1,VxVector a2,VxVector pt){
	VxVector l = a2 - a1;
	l.Normalize();
	float length = VxVectorInnerProduct(l,pt - a1);
	VxVector mid = a1 + length * (a2 - a1);
	return (pt - mid).Magnitude();
}

static bool SegmentIntersection(VxVector a1,
    VxVector a2,VxVector b1,VxVector b2,VxVector ret){
	VxVector a = a2 - a1;
	VxVector lineb1 = b1 - a1;
	VxVector lineb2 = b2 - a1;
	float sign = VxVectorInnerProduct(VxVectorCrossProduct(a,lineb1),
		VxVectorCrossProduct(a,lineb2));
	if(sign >= 0)
		return false;

	float d1 = distance2line(a1,a2,b1);
	float d2 = distance2line(a1,a2,b2);

	float sum = d1 + d2;

	ret = d1 / sum * a1 + d2 / sum * a2;
	return true;
}

static void ShowVxVector(CKContext* context,const VxVector& v){
	//context->OutputToConsoleEx("x: %f,y: %f,z %f",v.x,v.y,v.z);
}

// in range [lb,ub]
static bool InTheRange(float lb,float ub,float value){
	if(lb > ub)
		throw std::string("lower bound must be greater than upper bound!");
	if(value < lb)
		return false;
	else if(value > ub)
		return false;
	else
		return true;
}

#define DEBUGBREAK {context->OutputToConsoleEx("%s , %d",__FILE__,__LINE__);throw std::string("break");}

//bool VxVevtorIntersect

#endif