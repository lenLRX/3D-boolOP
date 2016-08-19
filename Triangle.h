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

enum InclusionRelation{
	Out = 0,
	In = 1,
	OnTheFace = 2
};

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

	bool VxVevtorIntersectTriangle(VxVector pt,VxVector line,bool& onThePlane,bool& isParallel,VxVector& IntersectionPoint = VxVector(1,1,1)){
		bool _onThePlane = false;
		bool _parallel = false;
		if(plane.RayIntersectTest(pt,line,_onThePlane,_parallel,IntersectionPoint)){
			VxVector AP = IntersectionPoint - v[0];
			float gamma = VxVectorInnerProduct(AP, _v);
			if(gamma >= 0 - 0.001 && gamma <= 1 + 0.001) {
				float beta = VxVectorInnerProduct(AP, _w);
				if(beta >= 0 - 0.001  && beta <= 1 + 0.001) {
                    float alpha = 1 - gamma - beta;
					if(alpha >= 0 - 0.001  && alpha <= 1 + 0.001){
						onThePlane = _onThePlane;
						isParallel = _parallel;
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

class PointInTriangle
{
public:
	PointInTriangle(Triangle T,VxVector pt){
		point = pt;
		pointInTriangle = T.PositionOfPointInTriangle(pt);
		if(fabs(pointInTriangle.v[0]) < 0.001)//alpha == 0
		{
			OnTheEdge = 0;
		}
		else if(fabs(pointInTriangle.v[1]) < 0.001)
		{
			OnTheEdge = 1;
		}
		else if(fabs(pointInTriangle.v[2]) < 0.001)
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

struct TriangleIntersection{
	TriangleIntersection(){}
	TriangleIntersection(const TriangleIntersection& TI){
		T1 = TI.T1;
		T2 = TI.T2;
		V1 = TI.V1;
		V2 = TI.V2;
		T1valid = TI.T1valid;
		T2valid = TI.T2valid;
	}
	TriangleIntersection(Triangle _T1,Triangle _T2,VxVector _V1,VxVector _V2,CKContext* context = NULL):
        T1(_T1),T2(_T2),V1(_V1),V2(_V2),T1valid(true),T2valid(true){

			if(context){
			PointInTriangle p1(T1,V1);
					if(!p1.valid()){
						context->OutputToConsoleEx("x: %f,y: %f,z %f",p1.pointInTriangle.x,p1.pointInTriangle.y,p1.pointInTriangle.z);
						DEBUGBREAK
					}

					PointInTriangle p2(T1,V2);
					if(!p2.valid()){
						context->OutputToConsoleEx("x: %f,y: %f,z %f",p2.pointInTriangle.x,p2.pointInTriangle.y,p2.pointInTriangle.z);
						DEBUGBREAK
					}

					PointInTriangle p3(T2,V1);
					if(!p3.valid()){
						context->OutputToConsoleEx("x: %f,y: %f,z %f",p3.pointInTriangle.x,p3.pointInTriangle.y,p3.pointInTriangle.z);
						DEBUGBREAK
					}

					PointInTriangle p4(T2,V2);
					if(!p4.valid()){
						context->OutputToConsoleEx("x: %f,y: %f,z %f",p4.pointInTriangle.x,p4.pointInTriangle.y,p4.pointInTriangle.z);
						DEBUGBREAK
					}
		}

		}
	Triangle T1,T2;
	VxVector V1,V2;
	bool T1valid;
	bool T2valid;
};




typedef std::vector<TriangleIntersection> Intersections;

static Intersections IntersectInplane(CKContext* context,Triangle T1,Triangle T2){
	//context->OutputToConsoleEx("IntersectInplane");
	Intersections ret;
	bool T1inT2[3] = {false};
	bool T2inT1[3] = {false};

	int T1inT2Count = 0;
    int T2inT1Count = 0;
	for(int i = 0;i < 3;i++){
		PointInTriangle p1(T2,T1.v[i]);
		PointInTriangle p2(T1,T2.v[i]);
		if(p1.valid() && p1.OnTheEdge < 0){
			T1inT2[i] = true;
			T1inT2Count++;
		}

		if(p2.valid() && p2.OnTheEdge < 0){
			T2inT1[i] = true;
			T2inT1Count++;
		}
	}

	if(3 == T1inT2Count){
		if(T2inT1Count){
			DEBUGBREAK
			return Intersections();
		}else{
			for(int i = 0;i < 3;i++){
				TriangleIntersection intersection(T1,T2,T1.v[i],T1.v[(i + 1) % 3],context);
				intersection.T1valid = false;
				intersection.T2valid = true;
				ret.push_back(intersection);
			}
		}
	}else if(3 == T2inT1Count){
		if(T1inT2Count){
			DEBUGBREAK
			return Intersections();
		}else{
			for(int i = 0;i < 3;i++){
				TriangleIntersection intersection(T1,T2,T2.v[i],T2.v[(i + 1) % 3],context);
				intersection.T1valid = false;
				intersection.T2valid = false;
				ret.push_back(intersection);
			}
		}
	}else if(0 == T1inT2Count){
		if(0 == T2inT1Count){
			/*
			for(int i = 0;i < 3;i++){
				std::vector<VxVector> x;
				for(int j = 0; j < 3;j++){
					VxVector cross;
					bool b = SegmentIntersection(T1.v[i],T1.v[(i + 1) % 3],T2.v[j],T2.v[(j + 1) % 3],cross);
					if(b){
						x.push_back(cross);
					}
				}

				if(x.size() >= 2){
					TriangleIntersection intersection(T1,T2,x[0],x[1],context);
					intersection.T1valid = false;
					intersection.T2valid = true;
					ret.push_back(intersection);
					return ret;
				}
			}
			return ret;
			*/
			return Intersections();
		}else if(1 == T2inT1Count){
			int theT2PointinT1 = -1;//let it crash
			for(int i = 0;i < 3;i++){
				if(T2inT1[i]){
					theT2PointinT1 = i;
					break;
				}
			}

			const VxVector& T2origin = T2.v[theT2PointinT1];
			const VxVector& T2e1 = T2.v[(theT2PointinT1 + 1) % 3];
			const VxVector& T2e2 = T2.v[(theT2PointinT1 + 2) % 3];
			bool T2gote1 = false;
			bool T2gote2 = false;

			VxVector cross1;
			VxVector cross2;

			for(int i = 0;i < 3;i++){
				if(!T2gote1){
					bool b = SegmentIntersection(T2origin,T2e1,T1.v[i],T1.v[(i + 1) % 3],cross1);
					if(b){
						T2gote1 = true;
					}
				}
				if(!T2gote2){
					bool b = SegmentIntersection(T2origin,T2e2,T1.v[i],T1.v[(i + 1) % 3],cross2);
					if(b){
						T2gote2 = true;
					}
				}
			}//for

			if(!(T2gote1 && T2gote2))
				return ret;
			
			TriangleIntersection i1(T1,T2,cross1,cross2,context);// in T2
			i1.T1valid = false;
			i1.T2valid = true;
			

			TriangleIntersection i2(T1,T2,cross1,T2origin,context);
			i2.T1valid = false;
			i2.T2valid = false;

			TriangleIntersection i3(T1,T2,cross2,T2origin,context);
			i3.T1valid = false;
			i3.T2valid = false;

			ret.push_back(i1);
			ret.push_back(i2);
			ret.push_back(i3);
			return ret;
		}//if(1 == T2inT1Count)
		else if(2 == T2inT1Count){
			int theT2PointNotinT1 = -1;
			for(int i = 0;i < 3;i++){
				if(!T2inT1[i]){
					theT2PointNotinT1 = i;
					break;
				}
			}

			const VxVector& T2origin = T2.v[theT2PointNotinT1];
			const VxVector& T2e1 = T2.v[(theT2PointNotinT1 + 1) % 3];
			const VxVector& T2e2 = T2.v[(theT2PointNotinT1 + 2) % 3];
			bool T2gote1 = false;
			bool T2gote2 = false;

			VxVector cross1;
			VxVector cross2;

			for(int i = 0;i < 3;i++){
				if(!T2gote1){
					bool b = SegmentIntersection(T2origin,T2e1,T1.v[i],T1.v[(i + 1) % 3],cross1);
					if(b){
						T2gote1 = true;
					}
				}
				if(!T2gote2){
					bool b = SegmentIntersection(T2origin,T2e2,T1.v[i],T1.v[(i + 1) % 3],cross2);
					if(b){
						T2gote2 = true;
					}
				}
			}//for

			if(!(T2gote1 && T2gote2))
				return ret;
			
			TriangleIntersection i1(T1,T2,cross1,cross2,context);
			i1.T1valid = false;
			i1.T2valid = true;
			

			TriangleIntersection i2(T1,T2,T2e1,T2e2,context);
			i2.T1valid = false;
			i2.T2valid = false;

			TriangleIntersection i3(T1,T2,T2e1,cross1,context);
			i3.T1valid = false;
			i3.T2valid = false;

			TriangleIntersection i4(T1,T2,T2e2,cross2,context);
			i4.T1valid = false;
			i4.T2valid = false;

			ret.push_back(i1);
			ret.push_back(i2);
			ret.push_back(i3);
			ret.push_back(i4);
			return ret;
		}//if(2 == T2inT1Count)
		else
			DEBUGBREAK
	}else if(1 == T1inT2Count){
		int theT1PointinT2 = -1;//let it crash
		for(int i = 0;i < 3;i++){
			if(T1inT2[i]){
				theT1PointinT2 = i;
				break;
			}
		}

		const VxVector& origin = T1.v[theT1PointinT2];
		const VxVector& e1 = T1.v[(theT1PointinT2 + 1) % 3];
		const VxVector& e2 = T1.v[(theT1PointinT2 + 2) % 3];
		bool gote1 = false;
		bool gote2 = false;

		VxVector cross1;
		int cross1No = -1;
		VxVector cross2;
		int cross2No = -1;
		for(int i = 0;i < 3;i++){
			if(!gote1){
				bool b = SegmentIntersection(origin,e1,T2.v[i],T2.v[(i + 1) % 3],cross1);
				if(T2inT1[i]){
					cross1No = i + 1;
				}else
				    cross1No = i;
				if(b){
					gote1 = true;
				}
			}
			if(!gote2){
				bool b = SegmentIntersection(origin,e2,T2.v[i],T2.v[(i + 1) % 3],cross2);
				if(T2inT1[i]){
					cross2No = i + 1;
				}else{
					cross2No = i;
				}
				if(b){
					gote2 = true;
				}
			}
		}//for

		if(!(gote1 && gote2))
				return ret;

		if(0 == T2inT1Count){

			TriangleIntersection i1(T1,T2,cross1,cross2,context);
			i1.T1valid = false;
			i1.T2valid = false;
			
			
			TriangleIntersection i2(T1,T2,cross1,origin,context);
			i2.T1valid = false;
			i2.T2valid = true;

			TriangleIntersection i3(T1,T2,cross2,origin,context);
			i3.T1valid = false;
			i3.T2valid = true;
			

			ret.push_back(i1);
			ret.push_back(i2);
			ret.push_back(i3);
			return ret;
		}//if(0 == T2inT1Count)
		else if(1 == T2inT1Count){
			int theT2PointinT1;
			for(int i = 0;i < 3;i++){
				if(T2inT1[i]){
					theT2PointinT1 = i;
					break;
				}
			}

			const VxVector& T2origin = T2.v[theT2PointinT1];
			TriangleIntersection i1(T1,T2,T2origin,cross1,context);
			i1.T1valid = false;
			i1.T2valid = false;

			TriangleIntersection i2(T1,T2,T2origin,cross2,context);
			i2.T1valid = false;
			i2.T2valid = false;

			
			TriangleIntersection i3(T1,T2,origin,cross1,context);
			i3.T1valid = false;
			i3.T2valid = true;

			TriangleIntersection i4(T1,T2,origin,cross2,context);
			i4.T1valid = false;
			i4.T2valid = true;
			

			ret.push_back(i1);
			ret.push_back(i2);
			ret.push_back(i3);
			ret.push_back(i4);
			return ret;
		}//if(1 == T2inT1Count)
		else if(2 == T2inT1Count){
			int theT2PointNotinT1 = -1;
			for(int i = 0;i < 3;i++){
				if(!T2inT1[i]){
					theT2PointNotinT1 = i;
					break;
				}
			}

			const VxVector& T2origin = T2.v[theT2PointNotinT1];

			
			TriangleIntersection i1(T1,T2,T2origin,cross1,context);
			i1.T1valid = false;
			i1.T2valid = true;

			TriangleIntersection i2(T1,T2,T2origin,cross2,context);
			i2.T1valid = false;
			i2.T2valid = true;
			

			TriangleIntersection i3(T1,T2,T2.v[(theT2PointNotinT1 + 1) % 3],T2.v[(theT2PointNotinT1 + 2) % 3],context);
			i3.T1valid = false;
			i3.T2valid = false;


			TriangleIntersection i4;
			i4.T1 = T1;
			i4.T2 = T2;
			i4.T1valid = false;
			i4.T2valid = false;

			TriangleIntersection i5;
			i5.T1 = T1;
			i5.T2 = T2;
			i5.T1valid = false;
			i5.T2valid = false;

			i4.V1 = cross1;
			i4.V2 = T2.v[(theT2PointNotinT1 + 1) % 3];
			i5.V1 = cross2;
			i5.V2 = T2.v[(theT2PointNotinT1 + 2) % 3];
			if(cross1No != (theT2PointNotinT1 + 1) % 3){
				std::swap(i4.V2,i5.V2);
			}

			ret.push_back(i1);
			ret.push_back(i2);
			ret.push_back(i3);
            //ret.push_back(i4);
			//ret.push_back(i5);
			return ret;
		}//if(2 == T2inT1Count)
		else
			DEBUGBREAK
	}//1 == T1inT2Count
	else if(2 == T1inT2Count){
		if(0 == T2inT1Count){
			int theT1PointNotinT2 = -1;
			for(int i = 0;i < 3;i++){
				if(!T1inT2[i]){
					theT1PointNotinT2 = i;
					break;
				}
			}

			const VxVector& T1origin = T1.v[theT1PointNotinT2];
			const VxVector& T1e1 = T1.v[(theT1PointNotinT2 + 1) % 3];
			const VxVector& T1e2 = T1.v[(theT1PointNotinT2 + 2) % 3];
			bool T1gote1 = false;
			bool T1gote2 = false;

			VxVector cross1;
			VxVector cross2;

			for(int i = 0;i < 3;i++){
				if(!T1gote1){
					bool b = SegmentIntersection(T1origin,T1e1,T2.v[i],T2.v[(i + 1) % 3],cross1);
					if(b){
						T1gote1 = true;
					}
				}
				if(!T1gote2){
					bool b = SegmentIntersection(T1origin,T1e2,T2.v[i],T2.v[(i + 1) % 3],cross2);
					if(b){
						T1gote2 = true;
					}
				}
			}//for

			if(!(T1gote1 && T1gote2))
				return ret;

			TriangleIntersection i1(T1,T2,cross1,cross2,context);
			i1.T1valid = false;
			i1.T2valid = false;

			
			TriangleIntersection i2(T1,T2,cross1,T1e1,context);
			i2.T1valid = false;
			i2.T2valid = true;

			TriangleIntersection i3(T1,T2,cross2,T1e2,context);
			i3.T1valid = false;
			i3.T2valid = true;

			TriangleIntersection i4(T1,T2,T1e1,T1e2,context);
			i4.T1valid = false;
			i4.T2valid = true;
			

			ret.push_back(i1);
			ret.push_back(i2);
			ret.push_back(i3);
			ret.push_back(i4);
			return ret;

		}//if(0 == T2inT1Count)
		else if(1 == T2inT1Count){
			int theT1PointNotinT2 = -1;
			for(int i = 0;i < 3;i++){
				if(!T1inT2[i]){
					theT1PointNotinT2 = i;
					break;
				}
			}

			int theT2PointinT1;
			for(int i = 0;i < 3;i++){
				if(T2inT1[i]){
					theT2PointinT1 = i;
					break;
				}
			}

			const VxVector& T1origin = T1.v[theT1PointNotinT2];
			const VxVector& T1e1 = T1.v[(theT1PointNotinT2 + 1) % 3];
			const VxVector& T1e2 = T1.v[(theT1PointNotinT2 + 2) % 3];
			bool T1gote1 = false;
			bool T1gote2 = false;

			VxVector cross1;
			VxVector cross2;

			for(int i = 0;i < 3;i++){
				if(!T1gote1){
					bool b = SegmentIntersection(T1origin,T1e1,T2.v[i],T2.v[(i + 1) % 3],cross1);
					if(b){
						T1gote1 = true;
					}
				}
				if(!T1gote2){
					bool b = SegmentIntersection(T1origin,T1e2,T2.v[i],T2.v[(i + 1) % 3],cross2);
					if(b){
						T1gote2 = true;
					}
				}
			}//for

			if(!(T1gote1 && T1gote2))
				return ret;
			
			TriangleIntersection i1(T1,T2,cross1,cross2,context);
			i1.T1valid = false;
			i1.T2valid = true;

			TriangleIntersection i2(T1,T2,cross1,T1e1,context);
			i2.T1valid = false;
			i2.T2valid = true;

			TriangleIntersection i3(T1,T2,cross2,T1e2,context);
			i3.T1valid = false;
			i3.T2valid = true;
			

			TriangleIntersection i4(T1,T2,cross2,T2.v[theT2PointinT1],context);
			i4.T1valid = false;
			i4.T2valid = false;

			TriangleIntersection i5(T1,T2,cross1,T2.v[theT2PointinT1],context);
			i5.T1valid = false;
			i5.T2valid = false;

			ret.push_back(i1);
			ret.push_back(i2);
			ret.push_back(i3);
			ret.push_back(i4);
			ret.push_back(i5);
			return ret;
		}//if(1 == T2inT1Count)
		else
			DEBUGBREAK
	}//if(2 == T1inT2Count)
	else
		DEBUGBREAK

}

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
				if(line.SquareMagnitude() > 0.01f){
					VxVector x = VxVectorInnerProduct(line,n1);
					//check if line is on the plane
					if(x.SquareMagnitude() < 0.01f){
						OnSamePlane = true;
					}
					break;
				}
			}

			if(!OnSamePlane){
				return Intersections();
			}else{
				return IntersectInplane(context,T1,T2);
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


#endif