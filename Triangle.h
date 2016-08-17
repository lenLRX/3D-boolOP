#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include <algorithm>
#include <vector>
#include "CKAll.h"

#include "Plane.h"

static char* errorMsg = "Hit Edge!";

//       0
// 2           1
//
//1                2
//
//        0

//local function

static void __ShowVxVector(CKContext* context,const VxVector& v){
	//context->OutputToConsoleEx("x: %f,y: %f,z %f",v.x,v.y,v.z);
}

class Triangle
{
public:
	Triangle(){}

	Triangle(VxVector _v1,VxVector _v2,VxVector _v3,CKVINDEX _faceIndex):
	faceIndex(_faceIndex){
		v[0] = _v1;
		v[1] = _v2;
		v[2] = _v3;
		VxVector V01 = v[1] - v[0]; //b
		VxVector V02 = v[2] - v[0]; //a
		norm = VxVectorCrossProduct(V02,V01);
		norm.Normalize();
		plane = Plane(norm,v[0]);

		//_v = Vector3D.GetCrossProduct(this.Plane.N, this.B - this.A);
		_v = VxVectorCrossProduct(plane.norm,V01);
		//_v /= Vector3D.GetDotProduct(this.C - this.A, _v);
		_v = _v / VxVectorInnerProduct(V02,_v);

        //_w = Vector3D.GetCrossProduct(this.Plane.N, this.C - this.A);
        _w = VxVectorCrossProduct(plane.norm,V02);
        //_w /= Vector3D.GetDotProduct(this.B - this.A, _w);
        _w = _w / VxVectorInnerProduct(V01,_w);
	}



	//@param pt point that determines line
    //@param line direction of line

	bool VxVevtorIntersectTriangle(VxVector pt,VxVector line,VxVector& IntersectionPoint = VxVector(1,1,1)){
		if(plane.RayIntersectTest(pt,line,IntersectionPoint)){
			VxVector AP = IntersectionPoint - v[0];
			float gamma = VxVectorInnerProduct(AP, _v);
			if(gamma >= 0 - 0.001 && gamma <= 1 + 0.001) {
				float beta = VxVectorInnerProduct(AP, _w);
				if(beta >= 0 - 0.001 && beta <= 1 + 0.001) {
                    float alpha = 1 - gamma - beta;
					if(alpha >= 0 - 0.001 && alpha <=1 + 0.001){
						return true;
					}
                }
			}
		}
		return false;

	}

	bool VxVevtorIntersectTriangle(VxVector pt,VxVector line,int flag,VxVector& IntersectionPoint = VxVector(1,1,1)){
		if(plane.RayIntersectTest(pt,line,flag,IntersectionPoint)){
			VxVector AP = IntersectionPoint - v[0];
			float gamma = VxVectorInnerProduct(AP, _v);
			if(gamma >= 0 - 0.001 && gamma <= 1 + 0.001) {
				float beta = VxVectorInnerProduct(AP, _w);
				if(beta >= 0 - 0.001 && beta <= 1 + 0.001) {
                    float alpha = 1 - gamma - beta;
					if(alpha >= 0 - 0.001 && alpha <=1 + 0.001){
						return true;
					}
                }
			}
		}
		return false;

	}

	VxVector PositionOfPointInTriangle(VxVector pt){
		float gamma = VxVectorInnerProduct(pt - v[0], _v);
		float beta = VxVectorInnerProduct(pt - v[0], _w);
		float alpha = 1 - gamma - beta;

		return VxVector(alpha,beta,gamma);
	}

	VxVector norm;
	bool contains();

	VxVector v[3];
	bool visible[3];

	CKVINDEX faceIndex;
	Plane plane;
private:
	//precomputations : http://www.tuicool.com/articles/rUrMJvi
	VxVector _v;
	VxVector _w;
};

struct TriangleIntersection{
	TriangleIntersection(Triangle _T1,Triangle _T2,VxVector _V1,VxVector _V2):
        T1(_T1),T2(_T2),V1(_V1),V2(_V2),T1valid(true),T2valid(true){;}
	Triangle T1,T2;
	VxVector V1,V2;
	bool T1valid;
	bool T2valid;
};


class PointInTriangle
{
public:
	PointInTriangle(Triangle T,VxVector pt){
		point = pt;
		pointInTriangle = T.PositionOfPointInTriangle(pt);
		if(AlmostEqualZero(pointInTriangle.v[0]))//alpha == 0
		{
			OnTheEdge = 0;
		}
		else if(AlmostEqualZero(pointInTriangle.v[1]))
		{
			OnTheEdge = 1;
		}
		else if(AlmostEqualZero(pointInTriangle.v[2]))
		{
			OnTheEdge = 2;
		}
		else
			OnTheEdge = -1;
	}
	PointInTriangle(){}

	bool valid(){
		for(int i = 0;i < 3;i++){
			if(pointInTriangle.v[i] > 1.0f + 0.001 || pointInTriangle.v[i] < 0 - 0.001)
				return false;
		}
		return true;
	}

	VxVector point;
	VxVector pointInTriangle;//¦ÁA + ¦ÂB + ¦ÃC A v[0] B v[1] C v[2]
	int OnTheEdge;//-1: not on the edge,0:On BC 1:On AC 2:On AB
};

typedef std::vector<TriangleIntersection> Intersections;

static Intersections IntersectWith(CKContext* context,Triangle T1,Triangle T2){
		VxVector resultPt1;
		VxVector resultPt2;

		Intersections ret;

		//All vertices of this triangle is on one side of T2
		VxVector N2 = VxVectorCrossProduct(T2.v[1] - T2.v[0],T2.v[2] - T2.v[0]);
		VxVector N1 = VxVectorCrossProduct(T1.v[1] - T1.v[0],T1.v[2] - T1.v[0]);

		VxVector n1 = N1;
		n1.Normalize();
		VxVector n2 = N2;
		n2.Normalize();


		//check if two triangle parallels
		if(SameVertex(n1,n2) || SameVertex(n1,-n2)){
			bool OnSamePlane = false;
			//check if they are on the same plane
			for(int i = 0; i < 3;i++){
				VxVector line = T1.v[0] - T2.v[i];
				//check if two point is too close
				if(line.SquareMagnitude() > 0.001f){
					VxVector x = VxVectorCrossProduct(line,n1);
					//check if line is on the plane
					if(x.SquareMagnitude() < 0.0001f){
						OnSamePlane = true;
					}
					break;
				}
			}

			if(!OnSamePlane){
				return Intersections();
			}else{
				return Intersections();
			}
		}



		float d2 = - VxVectorInnerProduct(N2,T2.v[0]);
		float d2v[3];
		for(int i = 0;i < 3;i++){
			d2v[i] = VxVectorInnerProduct(N2,T1.v[i]) + d2;
		}

		bool d2AllNeg = d2v[0] < 0 && d2v[1] <0 && d2v[2] < 0;
		bool d2AllPos = d2v[0] > 0 && d2v[1] >0 && d2v[2] > 0;

		if(d2AllNeg || d2AllPos)
			return Intersections();

		float d1 = - VxVectorInnerProduct(N1,T1.v[0]);
		float d1v[3];
		for(int i = 0;i < 3;i++){
			d1v[i] = VxVectorInnerProduct(N1,T2.v[i]) + d1;
		}

		bool d1AllNeg = d1v[0] < 0 && d1v[1] <0 && d1v[2] < 0;
		bool d1AllPos = d1v[0] > 0 && d1v[1] >0 && d1v[2] > 0;

		if(d1AllNeg || d1AllPos)
			return Intersections();

		VxVector D = VxVectorCrossProduct(N1,N2);

		int T1indexTrans[3];//T1index tranformation

		if(d2v[0] * d2v[1] > 0)//T1 v0,v1 on one side
		{
			T1indexTrans[0] = 0;
			T1indexTrans[1] = 2;
			T1indexTrans[2] = 1;
			//context->OutputToConsoleEx("T1 v0,v1 on one side");
		}
		else if(d2v[0] * d2v[2] > 0)//T1 v0,v2 on one side
		{
			T1indexTrans[0] = 0;
			T1indexTrans[1] = 1;
			T1indexTrans[2] = 2;
			//context->OutputToConsoleEx("T1 v0,v2 on one side");
		}
		else//T1 v1,v2 on one side
		{
			T1indexTrans[0] = 1;
			T1indexTrans[1] = 0;
			T1indexTrans[2] = 2;
			//context->OutputToConsoleEx("T1 v1,v2 on one side");
		}

		float pvT1[3];

		for(int i = 0;i < 3;i++){
			pvT1[T1indexTrans[i]] = VxVectorInnerProduct(D,T1.v[T1indexTrans[i]]);
		}

		//T1

		float t1 = pvT1[T1indexTrans[0]] + 
			(pvT1[T1indexTrans[1]] - pvT1[T1indexTrans[0]]) * d2v[T1indexTrans[0]]/(d2v[T1indexTrans[0]] - d2v[T1indexTrans[1]]);

		float t2 = pvT1[T1indexTrans[2]] + 
			(pvT1[T1indexTrans[1]] - pvT1[T1indexTrans[2]]) * d2v[T1indexTrans[2]]/(d2v[T1indexTrans[2]] - d2v[T1indexTrans[1]]);

		if(t1 == t2)
			return Intersections();

		if(t1 > t2)
		{
			std::swap(t1,t2);
			std::swap(T1indexTrans[0],T1indexTrans[2]);
			//context->OutputToConsoleEx("swap 12");
		}
			

		int T2indexTrans[3];//T1index tranformation

		if(d1v[0] * d1v[1] > 0)//T2 v0,v1 on one side
		{
			T2indexTrans[0] = 0;
			T2indexTrans[1] = 2;
			T2indexTrans[2] = 1;
			//context->OutputToConsoleEx("T2 v0,v1 on one side");
		}
		else if(d1v[0] * d1v[2] > 0)//T2 v0,v2 on one side
		{
			T2indexTrans[0] = 0;
			T2indexTrans[1] = 1;
			T2indexTrans[2] = 2;
			//context->OutputToConsoleEx("T2 v0,v2 on one side");
		}
		else//T2 v1,v2 on one side
		{
			T2indexTrans[0] = 1;
			T2indexTrans[1] = 0;
			T2indexTrans[2] = 2;
			//context->OutputToConsoleEx("T2 v1,v2 on one side");
		}

		float pvT2[3];

		for(int i = 0;i < 3;i++){
			pvT2[T2indexTrans[i]] = VxVectorInnerProduct(D,T2.v[T2indexTrans[i]]);
		}

		float t3 = pvT2[T2indexTrans[0]] + 
			(pvT2[T2indexTrans[1]] - pvT2[T2indexTrans[0]]) * d1v[T2indexTrans[0]]/(d1v[T2indexTrans[0]] - d1v[T2indexTrans[1]]);

		float t4 = pvT2[T2indexTrans[2]] + 
			(pvT2[T2indexTrans[1]] - pvT2[T2indexTrans[2]]) * d1v[T2indexTrans[2]]/(d1v[T2indexTrans[2]] - d1v[T2indexTrans[1]]);

		if(t3 == t4)
			return Intersections();

		if(t3 > t4)
		{
			std::swap(t3,t4);
			std::swap(T2indexTrans[0],T2indexTrans[2]);
			//context->OutputToConsoleEx("swap 34");
		}
			

		if(t1 >= t4)// t3----t4....t1----t2 reject!
			return Intersections();

		if(t3 >= t2)// t1----t2....t3----t4 reject!
			return Intersections();

		if(t1 < t3){//t1----t3
			//t2<t3:we had proved above
			if(t2 > t4){//t1---t3--t4---t2
				//so the overlap part is t3---t4
				if(!T1.VxVevtorIntersectTriangle(T2.v[T2indexTrans[0]],
					T2.v[T2indexTrans[1]] - T2.v[T2indexTrans[0]],resultPt1))
					return Intersections();
				if(!T1.VxVevtorIntersectTriangle(T2.v[T2indexTrans[2]],
					T2.v[T2indexTrans[1]] - T2.v[T2indexTrans[2]],resultPt2))
					return Intersections();
				if(SameVertex(resultPt1,resultPt2))
					return Intersections();
				//context->OutputToConsoleEx("type 1");
				ret.push_back(TriangleIntersection(T1,T2,resultPt1,resultPt2));
				return ret;
			}
			else{//t1---t3---t2--t4
				//so the overlap part is t3---t2
				if(!T1.VxVevtorIntersectTriangle(T2.v[T2indexTrans[0]],
					T2.v[T2indexTrans[1]] - T2.v[T2indexTrans[0]],resultPt1))
					return Intersections();
				if(!T2.VxVevtorIntersectTriangle(T1.v[T1indexTrans[2]],
					T1.v[T1indexTrans[1]] - T1.v[T1indexTrans[2]],resultPt2))
					return Intersections();
				if(SameVertex(resultPt1,resultPt2))
					return Intersections();
				//context->OutputToConsoleEx("type 2");
				ret.push_back(TriangleIntersection(T1,T2,resultPt1,resultPt2));
				return ret;
			}
		}else{//t3----t1
			//t1<t4:we had proved above
			if(t2>t4){//t3--t1---t4--t2
				//so the overlap part is t1---t4
				if(!T2.VxVevtorIntersectTriangle(T1.v[T1indexTrans[0]],
					T1.v[T1indexTrans[1]] - T1.v[T1indexTrans[0]],resultPt1))
					return Intersections();
				if(!T1.VxVevtorIntersectTriangle(T2.v[T2indexTrans[2]],
					T2.v[T2indexTrans[1]] - T2.v[T2indexTrans[2]],resultPt2))
					return Intersections();
				if(SameVertex(resultPt1,resultPt2))
					return Intersections();
				//context->OutputToConsoleEx("type 3");
				ret.push_back(TriangleIntersection(T1,T2,resultPt1,resultPt2));
				return ret;
			}else{//t3--t1---t2---t4
				//so the overlap part is t1---t2
				//T1V0 -> T2 -> T1V1
				if(!T2.VxVevtorIntersectTriangle(T1.v[T1indexTrans[0]],
					T1.v[T1indexTrans[1]] - T1.v[T1indexTrans[0]],resultPt1)){
					return Intersections();
					__ShowVxVector(context,resultPt1);
					__ShowVxVector(context,T1.v[T1indexTrans[0]]);
					__ShowVxVector(context,T1.v[T1indexTrans[1]]);
					__ShowVxVector(context,T1.v[T1indexTrans[2]]);
					__ShowVxVector(context,T1.v[T1indexTrans[1]] - T1.v[T1indexTrans[0]]);
					__ShowVxVector(context,T2.norm);
					__ShowVxVector(context,T2.v[T2indexTrans[0]]);
					__ShowVxVector(context,T2.v[T2indexTrans[1]]);
					__ShowVxVector(context,T2.v[T2indexTrans[2]]);
					DEBUGBREAK
					return Intersections();
					//DEBUGBREAK
				}
				if(!T2.VxVevtorIntersectTriangle(T1.v[T1indexTrans[2]],
					T1.v[T1indexTrans[1]] - T1.v[T1indexTrans[2]],resultPt2))
					return Intersections();//DEBUGBREAK
				if(SameVertex(resultPt1,resultPt2))
					return Intersections();
				//context->OutputToConsoleEx("type 4");
				ret.push_back(TriangleIntersection(T1,T2,resultPt1,resultPt2));
				return ret;
			}
		}


		//float pv0T1 = VxVectorInnerProduct(D,this.v[0]);

	}

static Intersections IntersectInplane(CKContext* context,Triangle T2){
		
}


#endif