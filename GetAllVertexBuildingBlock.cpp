//////////////////////////////////////////////////////////////////////////////////////////////////////////
//		            GetAllVertexBuildingBlock
//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "CKAll.h"
#include <stdio.h>
#include <set>
#include <string>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include <time.h>
#include "MyUtility.h"
#include "Plane.h"
#include "Triangle.h"
#include "Polygon.h"

CKObjectDeclaration	*FillBehaviorGetAllVertexBuildingBlockDecl();
CKERROR CreateGetAllVertexBuildingBlockProto(CKBehaviorPrototype **);
int GetAllVertexBuildingBlock(const CKBehaviorContext& BehContext);
int GetAllVertexBuildingBlockCallBack(const CKBehaviorContext& BehContext);

CKObjectDeclaration	*FillBehaviorGetAllVertexBuildingBlockDecl()
{
	CKObjectDeclaration *od = CreateCKObjectDeclaration("bool cut");
	
	od->SetType(CKDLL_BEHAVIORPROTOTYPE);
	od->SetVersion(0x00000001);
	od->SetCreationFunction(CreateGetAllVertexBuildingBlockProto);
	od->SetDescription("Enter your description here");
	od->SetCategory("DLVR");
	od->SetGuid(CKGUID(0x4acca9a7, 0xb22da652));
	od->SetAuthorGuid(CKGUID(0x56495254,0x4f4f4c53));
	od->SetAuthorName("545976176@qq.com");
	od->SetCompatibleClassId(CKCID_BEOBJECT);

	return od;
}

CKERROR CreateGetAllVertexBuildingBlockProto(CKBehaviorPrototype** pproto)
{
	CKBehaviorPrototype *proto = CreateCKBehaviorPrototype("bool cut");
	if (!proto) {
		return CKERR_OUTOFMEMORY;
	}

//---     Inputs declaration
	proto->DeclareInput("In0");

//---     Outputs declaration
	proto->DeclareOutput("Out0");


	proto->DeclareInParameter("cutter", CKPGUID_3DENTITY );

//----	Local Parameters Declaration

//----	Settings Declaration
	proto->SetBehaviorFlags((CK_BEHAVIOR_FLAGS)(CKBEHAVIOR_TARGETABLE));
	proto->SetBehaviorCallbackFct(GetAllVertexBuildingBlockCallBack, CKCB_BEHAVIORATTACH|CKCB_BEHAVIORDETACH|CKCB_BEHAVIORDELETE|CKCB_BEHAVIOREDITED|CKCB_BEHAVIORSETTINGSEDITED|CKCB_BEHAVIORLOAD|CKCB_BEHAVIORPRESAVE|CKCB_BEHAVIORPOSTSAVE|CKCB_BEHAVIORRESUME|CKCB_BEHAVIORPAUSE|CKCB_BEHAVIORRESET|CKCB_BEHAVIORRESET|CKCB_BEHAVIORDEACTIVATESCRIPT|CKCB_BEHAVIORACTIVATESCRIPT|CKCB_BEHAVIORREADSTATE, NULL);
	proto->SetFunction(GetAllVertexBuildingBlock);

	*pproto = proto;
	return CK_OK;
}

InclusionRelation PointInBody(VxVector pt,VxVector scale,CKMesh* mesh,int flag,CKContext* context,VxVector direction = VxVector(0,-1,-1),bool review = false){

	CKDWORD vStride=0;

    const BYTE* vxv = (BYTE*)mesh->GetModifierVertices(&vStride);
	int vcount = mesh->GetModifierVertexCount();

	int faceCount = mesh->GetFaceCount();
	
	const CKVINDEX* faceIndex = mesh->GetFacesIndices();

	VxVector line = direction;

	line = VxVector(0,1,0);//(float)(rand() + 1)/(float)RAND_MAX * 0.05

	//context->OutputToConsoleEx("line: %f,%f,%f",line.x,line.y,line.z);

	int IntersectionCount = 0;

	//context->OutputToConsoleEx("face Count %d , vcount %d", faceCount, vcount);

	std::set<VxVector,VxVectorLess> IntersectionPointSet;

	for(int i = 0; i<faceCount; i++ , faceIndex+=3){
		CKVINDEX f1 = *(faceIndex);
		CKVINDEX f2 = *(faceIndex+1);
		CKVINDEX f3 = *(faceIndex+2);
		//context->OutputToConsoleEx("faces %d,%d,%d",f1,f2,f3);

		VxVector v1 = *((VxVector*)(vxv + f1*vStride ));
	    VxVector v2 = *((VxVector*)(vxv + f2*vStride ));
	    VxVector v3 = *((VxVector*)(vxv + f3*vStride ));

		for(int j = 0;j < 3;j++){
			v1.v[j] *= scale.v[j];
			v2.v[j] *= scale.v[j];
			v3.v[j] *= scale.v[j];
		}

		//context->OutputToConsoleEx("v1: %f,%f,%f",v1.x,v1.y,v1.z);
		//context->OutputToConsoleEx("v2: %f,%f,%f",v2.x,v2.y,v2.z);
		//context->OutputToConsoleEx("v3: %f,%f,%f",v3.x,v3.y,v3.z);



		Triangle triangle(v1,v2,v3,f1);

		VxVector intersectionPoint;

		bool onThePlane = false;
		bool isParallel = false;
		bool b = triangle.VxVevtorIntersectTriangle(pt,line,onThePlane,isParallel,intersectionPoint);

		if(onThePlane && b)
			return OnTheFace;

		if(isParallel && b){
			return OnTheFace;
		}
		/*
		if(b){
			//context->OutputToConsoleEx("intersectionPoint: %f,%f,%f",
			//	intersectionPoint.x,intersectionPoint.y,intersectionPoint.z);
			if(IntersectionPointSet.count(intersectionPoint)){
				b = false;
			}
			else{
				IntersectionPointSet.insert(intersectionPoint);
			}
		}
		*/
		
		if(b)
			IntersectionCount++;

		/*
		if(b){
			PointInTriangle pointInTriangle(triangle,pt);
			if(pointInTriangle.OnTheEdge < 0)
			    IntersectionCount += 2;
			else
				IntersectionCount++;
	    }
		*/
	}

	//IntersectionCount = IntersectionCount / 2;

	if(IntersectionCount % 2 == 1)
		return In;
	else
		return Out;
}

std::vector<InclusionRelation> Mesh1PointsInMesh2(VxBbox box,CKMesh* mesh1,CKMesh* mesh2,VxVector displacement,VxVector scale1,VxVector scale2,int flag,CKContext* context){
	CKDWORD Mesh1vStride=0;
	CKDWORD Mesh2vStride=0;
	const BYTE* Mesh1Vertices = (BYTE*)mesh1->GetModifierVertices(&Mesh1vStride);
	const BYTE* Mesh2Vertices = (BYTE*)mesh2->GetModifierVertices(&Mesh2vStride);
	int Mesh1VCount = mesh1->GetModifierVertexCount();
	int Mesh2VCount = mesh2->GetModifierVertexCount();

	std::vector<InclusionRelation> ret;
	int count = 0;
	
	box = mesh2->GetLocalBox();

	for(int i = 0;i < 6;i++){
		box.v[i] *= 1.2f;
	}

	//context->OutputToConsoleEx("Max: %f %f %f",box.Max.x,box.Max.y,box.Max.z);
	//context->OutputToConsoleEx("Min: %f %f %f",box.Min.x,box.Min.y,box.Min.z);
	
	for(int i=0; i<Mesh1VCount ;i++,Mesh1Vertices+=Mesh1vStride) {

		VxVector TransformedPos = *(VxVector*)Mesh1Vertices;
		/*
		for(int j = 0;j < 3;j++){
			TransformedPos.v[j] *= scale1.v[j];
		}
		*/

		
		TransformedPos += displacement;
			

		int Count[3] = {0};

		for(int j = 0 ;j < 1;j++){
			InclusionRelation ir = PointInBody(TransformedPos,scale2,mesh2,flag,context);
			Count[ir]++;
		}

		int max = -1;
		int pos = -1;

		for(int j = 0 ;j < 3;j++){
			if(max < Count[j]){
				max = Count[j];
				pos = j;
			}
		}
		ret.push_back((InclusionRelation)pos);

		/*
		if(box.VectorIn(TransformedPos))
		    ret.push_back(PointInBody(TransformedPos,scale2,mesh2,flag,context));
		else
			ret.push_back(Out);
			*/
	}
	

	return ret;
}

std::vector<TriangleIntersection> IntersectedTrianglesOf2Mesh(VxBbox box,CKMesh* mesh1,CKMesh* mesh2,
															  VxVector displacement,VxVector scale1,VxVector scale2,CKContext* context){
	std::vector<TriangleIntersection> results;
	//context->OutputToConsoleEx("IntersectedTrianglesOf2Mesh");
	CKDWORD Mesh1vStride=0;
	CKDWORD Mesh2vStride=0;
	const BYTE* Mesh1Vertices = (BYTE*)mesh1->GetModifierVertices(&Mesh1vStride);
	const BYTE* Mesh2Vertices = (BYTE*)mesh2->GetModifierVertices(&Mesh2vStride);
	int Mesh1VCount = mesh1->GetModifierVertexCount();
	int Mesh2VCount = mesh2->GetModifierVertexCount();
	int Mesh1FaceCount = mesh1->GetFaceCount();
	int Mesh2FaceCount = mesh2->GetFaceCount();
	const CKVINDEX* faceIndex1 = mesh1->GetFacesIndices();
	const CKVINDEX* faceIndex2 = mesh2->GetFacesIndices();

	VxBbox box1 = mesh1->GetLocalBox();
	VxBbox box2 = mesh2->GetLocalBox();

	for(int i = 0;i < 6;i++){
		box1.v[i] *= 1.2f;
		box2.v[i] *= 1.2f;
	}


	int IntersectionCount = 0;

	//context->OutputToConsoleEx("mesh1 face %d mesh2 face %d",Mesh1FaceCount, Mesh2FaceCount);

	for(int i = 0;i < Mesh1FaceCount;i++){
		CKVINDEX T1f1 = *(faceIndex1 + 3 * i);
		CKVINDEX T1f2 = *(faceIndex1 + 3 * i + 1);
		CKVINDEX T1f3 = *(faceIndex1 + 3 * i + 2);

		//context->OutputToConsoleEx("mesh1 faces %d,%d,%d",T1f1,T1f2,T1f3);

		VxVector T1v1 = *((VxVector*)(Mesh1Vertices + T1f1*Mesh1vStride ));
	    VxVector T1v2 = *((VxVector*)(Mesh1Vertices + T1f2*Mesh1vStride ));
	    VxVector T1v3 = *((VxVector*)(Mesh1Vertices + T1f3*Mesh1vStride ));

		/*
		for(int k = 0;k < 3;k++){
			T1v1.v[k] *= scale1.v[k];
			T1v2.v[k] *= scale1.v[k];
			T1v3.v[k] *= scale1.v[k];
		}
		*/

		/*
		if(!(box2.VectorIn(T1v1 - displacement) || box2.VectorIn(T1v2  - displacement) || box2.VectorIn(T1v3  - displacement)))
			continue;
			*/

		Triangle T1(T1v1,T1v2,T1v3,i);

		for(int j = 0;j < Mesh2FaceCount;j++){
	        CKVINDEX T2f1 = *(faceIndex2 + 3 * j);
			CKVINDEX T2f2 = *(faceIndex2 + 3 * j + 1);
			CKVINDEX T2f3 = *(faceIndex2 + 3 * j + 2);

			//context->OutputToConsoleEx("mesh2 faces %d,%d,%d",T2f1,T2f2,T2f3);

			VxVector T2v1 = *((VxVector*)(Mesh2Vertices + T2f1*Mesh2vStride )) + displacement;
			VxVector T2v2 = *((VxVector*)(Mesh2Vertices + T2f2*Mesh2vStride )) + displacement;
			VxVector T2v3 = *((VxVector*)(Mesh2Vertices + T2f3*Mesh2vStride )) + displacement;

			/*
			if(!(box1.VectorIn(T2v1) || box1.VectorIn(T2v2) || box1.VectorIn(T2v3)))
			    continue;
				*/

			/*
			for(int k = 0;k < 3;k++){
			    T2v1.v[k] *= scale2.v[k];
			    T2v2.v[k] *= scale2.v[k];
			    T2v3.v[k] *= scale2.v[k];
		    }
			*/

			Triangle T2(T2v1,T2v2,T2v3,j);

			//TODO:check if it is in triangle

			Intersections intersections = IntersectWith(context,T1,T2);
			if(intersections.size()){
				IntersectionCount += intersections.size();
				
				for(size_t idx = 0;idx < intersections.size();idx++){
					/*
					PointInTriangle p1(T1,intersections[idx].V1);
					if(!p1.valid()){
						context->OutputToConsoleEx("x: %f,y: %f,z %f",p1.pointInTriangle.x,p1.pointInTriangle.y,p1.pointInTriangle.z);
						DEBUGBREAK
					}

					PointInTriangle p2(T1,intersections[idx].V2);
					if(!p2.valid()){
						context->OutputToConsoleEx("x: %f,y: %f,z %f",p2.pointInTriangle.x,p2.pointInTriangle.y,p2.pointInTriangle.z);
						DEBUGBREAK
					}

					PointInTriangle p3(T2,intersections[idx].V1);
					if(!p3.valid()){
						context->OutputToConsoleEx("x: %f,y: %f,z %f",p3.pointInTriangle.x,p3.pointInTriangle.y,p3.pointInTriangle.z);
						DEBUGBREAK
					}

					PointInTriangle p4(T2,intersections[idx].V2);
					if(!p4.valid()){
						context->OutputToConsoleEx("x: %f,y: %f,z %f",p4.pointInTriangle.x,p4.pointInTriangle.y,p4.pointInTriangle.z);
						DEBUGBREAK
					}
					*/

					results.push_back(intersections[idx]);
				}
			}
		}
	}
	return results;
}

void CutMesh1ByMesh2(CKMesh* mesh1,CKMesh* mesh2,VxBbox box1,VxBbox box2,
					 VxVector displacement,VxVector scale1,VxVector scale2,CKContext* context){
	context->OutputToConsoleEx("box1 Max: %f %f %f",box1.Max.x,box1.Max.y,box1.Max.z);
	context->OutputToConsoleEx("box1 Min: %f %f %f",box1.Min.x,box1.Min.y,box1.Min.z);
	context->OutputToConsoleEx("box2 Max: %f %f %f",box2.Max.x,box2.Max.y,box2.Max.z);
	context->OutputToConsoleEx("box2 Min: %f %f %f",box2.Min.x,box2.Min.y,box2.Min.z);
    box1.Intersect(box2);
    VxBbox box = box1;
    srand(time(NULL));
    std::vector<InclusionRelation> Mesh1InMesh2 = Mesh1PointsInMesh2(box,mesh1,mesh2,-displacement,scale1,scale2,1,context);
    std::vector<InclusionRelation> Mesh2InMesh1 = Mesh1PointsInMesh2(box,mesh2,mesh1, displacement,scale2,scale1,-1,context);

	/*
	size_t s = Mesh2InMesh1.second.size();
	for(size_t i = 0;i < s;i++){
		// visibility is inverse for Mesh2InMesh1
		Mesh2InMesh1.second[i] = !Mesh2InMesh1.second[i];
	}
	*/

	CKDWORD Mesh1vStride=0;
	CKDWORD Mesh2vStride=0;
	BYTE* Mesh1Vertices = (BYTE*)mesh1->GetModifierVertices(&Mesh1vStride);
	BYTE* Mesh2Vertices = (BYTE*)mesh2->GetModifierVertices(&Mesh2vStride);
	int Mesh1VCount = mesh1->GetModifierVertexCount();
	int Mesh2VCount = mesh2->GetModifierVertexCount();
	int Mesh1FaceCount = mesh1->GetFaceCount();
	int Mesh2FaceCount = mesh2->GetFaceCount();
	CKVINDEX* faceIndex1 = mesh1->GetFacesIndices();
	CKVINDEX* faceIndex2 = mesh2->GetFacesIndices();

	std::vector<TriangleIntersection> TriangleIntersections = IntersectedTrianglesOf2Mesh(box,mesh1,mesh2,displacement,scale1,scale2,context);

	size_t IntersectionSize = TriangleIntersections.size();

	std::vector<bool> Mesh1TriangleMarks(Mesh1FaceCount,false);
	std::vector<bool> Mesh2TriangleMarks(Mesh2FaceCount,false);

	std::vector<int> NewFaces;
	std::vector<VxVector> NewPoints;

	std::vector<Polygon> T1Polygons(Mesh1FaceCount);
	std::set<int> T1Triangles;
	std::vector<Polygon> T2Polygons(Mesh2FaceCount);
	std::set<int> T2Triangles;

	for(size_t i = 0; i < IntersectionSize ; i++){
		TriangleIntersection& TI = TriangleIntersections[i];

		Edge E;

		if(TI.T1valid){
			Mesh1TriangleMarks[TI.T1.faceIndex] = true;
			T1Polygons[TI.T1.faceIndex].triangle = TI.T1;

			T1Triangles.insert(TI.T1.faceIndex);

			for(int j = 0; j < 3;j++){
				switch(Mesh1InMesh2[*(faceIndex1 + 3 * TI.T1.faceIndex + j)]){
					case In:
						T1Polygons[TI.T1.faceIndex].triangle.visible[j] = false;
						break;
					case Out:
						T1Polygons[TI.T1.faceIndex].triangle.visible[j] = true;
						break;
					case OnTheFace:
						T1Polygons[TI.T1.faceIndex].triangle.visible[j] = false;
						/*
						InclusionRelation ir1 = PointInBody(T1Polygons[TI.T1.faceIndex].triangle.v[j],
							VxVector(1.0f,1.0f,1.0f),mesh2,1,context,
							T1Polygons[TI.T1.faceIndex].triangle.v[(j + 1) % 3 ] - T1Polygons[TI.T1.faceIndex].triangle.v[j],true);
						InclusionRelation ir2 = PointInBody(T1Polygons[TI.T1.faceIndex].triangle.v[j],
							VxVector(1.0f,1.0f,1.0f),mesh2,1,context,
							T1Polygons[TI.T1.faceIndex].triangle.v[(j + 2) % 3 ] - T1Polygons[TI.T1.faceIndex].triangle.v[j],true);
						if(OnTheFace == ir1 && OnTheFace == ir2){
							T1Polygons[TI.T1.faceIndex].triangle.visible[j] = false;
							//Mesh1InMesh2[*(faceIndex1 + 3 * TI.T1.faceIndex + j)] = Out;
						}
						if(In == ir1 || In == ir2){
							T1Polygons[TI.T1.faceIndex].triangle.visible[j] = false;
							//Mesh1InMesh2[*(faceIndex1 + 3 * TI.T1.faceIndex + j)] = In;
						}else{
							T1Polygons[TI.T1.faceIndex].triangle.visible[j] = false;
							//Mesh1InMesh2[*(faceIndex1 + 3 * TI.T1.faceIndex + j)] = Out;
						}
						*/
						break;
				}
				/*
				T1Polygons[TI.T1.faceIndex].triangle.visible[j]
				= !Mesh1InMesh2[*(faceIndex1 + 3 * TI.T1.faceIndex + j)];
				*/
			}

			int vi = T1Polygons[TI.T1.faceIndex].triangle.visible[0]
			+ T1Polygons[TI.T1.faceIndex].triangle.visible[1]
			+ T1Polygons[TI.T1.faceIndex].triangle.visible[3];

			if(vi == 0 || vi == 3){
				
				continue;
			}

			E.v1 = PointInTriangle(TI.T1,TI.V1);
			E.v2 = PointInTriangle(TI.T1,TI.V2);

			if(!E.v1.valid()){
				context->OutputToConsoleEx("x: %f,y: %f,z %f",E.v1.pointInTriangle.x,E.v1.pointInTriangle.y,E.v1.pointInTriangle.z);
				DEBUGBREAK
			}
			if(!E.v2.valid()){
				context->OutputToConsoleEx("x: %f,y: %f,z %f",E.v1.pointInTriangle.x,E.v1.pointInTriangle.y,E.v1.pointInTriangle.z);
				DEBUGBREAK
			}


			T1Polygons[TI.T1.faceIndex].Edges.push_back(E);
		}

		if(TI.T2valid){
			Mesh2TriangleMarks[TI.T2.faceIndex] = true;

			T2Polygons[TI.T2.faceIndex].triangle = TI.T2;

			T2Triangles.insert(TI.T2.faceIndex);

			for(int j = 0; j < 3;j++){
				switch(Mesh2InMesh1[*(faceIndex2 + 3 * TI.T2.faceIndex + j)]){
					case Out:
						T2Polygons[TI.T2.faceIndex].triangle.visible[j] = false;
						break;
					case In:
						T2Polygons[TI.T2.faceIndex].triangle.visible[j] = true;
						break;
					case OnTheFace:
						T2Polygons[TI.T2.faceIndex].triangle.visible[j] = false;
						/*
						InclusionRelation ir1 = PointInBody(T2Polygons[TI.T2.faceIndex].triangle.v[j],
							VxVector(1.0f,1.0f,1.0f),mesh1,1,context,
							T2Polygons[TI.T2.faceIndex].triangle.v[(j + 1) % 3 ] - T2Polygons[TI.T2.faceIndex].triangle.v[j],true);
						InclusionRelation ir2 = PointInBody(T2Polygons[TI.T2.faceIndex].triangle.v[j],
							VxVector(1.0f,1.0f,1.0f),mesh1,1,context,
							T2Polygons[TI.T2.faceIndex].triangle.v[(j + 2) % 3 ] - T2Polygons[TI.T2.faceIndex].triangle.v[j],true);
						if(OnTheFace == ir1 && OnTheFace == ir2){
							T2Polygons[TI.T2.faceIndex].triangle.visible[j] = false;
							//Mesh2InMesh1[*(faceIndex2 + 3 * TI.T2.faceIndex + j)] = Out;
						}
						else if(Out == ir1 || Out == ir2){
							T2Polygons[TI.T2.faceIndex].triangle.visible[j] = false;
							//Mesh2InMesh1[*(faceIndex2 + 3 * TI.T2.faceIndex + j)] = Out;
						}else{
							T2Polygons[TI.T2.faceIndex].triangle.visible[j] = false;
							//Mesh2InMesh1[*(faceIndex2 + 3 * TI.T2.faceIndex + j)] = In;
						}
						*/
						break;
				}
				/*
				T2Polygons[TI.T2.faceIndex].triangle.visible[j] 
				= !Mesh2InMesh1.second[*(faceIndex2 + 3 * TI.T2.faceIndex + j)];
				*/
			}

			int vi = T2Polygons[TI.T2.faceIndex].triangle.visible[0]
			+ T2Polygons[TI.T2.faceIndex].triangle.visible[1]
			+ T2Polygons[TI.T2.faceIndex].triangle.visible[2];

			if(vi == 0 || vi == 3){
				
				continue;
			}


			/*
			fprintf(fp,"visible: %d %d %d \n",T2Polygons[TI.T2.faceIndex].triangle.visible[0],
				T2Polygons[TI.T2.faceIndex].triangle.visible[1],T2Polygons[TI.T2.faceIndex].triangle.visible[2]);

			fprintf(fp,"IR: %d %d %d \n",Mesh2InMesh1[*(faceIndex2 + 3 * TI.T2.faceIndex + 0)],
					Mesh2InMesh1[*(faceIndex2 + 3 * TI.T2.faceIndex + 1)],Mesh2InMesh1[*(faceIndex2 + 3 * TI.T2.faceIndex + 2)]);

			for(int j = 0;j < 3;j++)
			    fprintf(fp,"T1.v[%d] : %f %f %f \n",j,TI.T1.v[j].x,TI.T1.v[j].y,TI.T1.v[j].z);

			for(int j = 0;j < 3;j++)
			    fprintf(fp,"T2.v[%d] : %f %f %f \n",j,TI.T2.v[j].x,TI.T2.v[j].y,TI.T2.v[j].z);

			fprintf(fp,"V1 x: %f,y: %f,z %f \n",TI.V1.x,TI.V1.y,TI.V1.z);
			fprintf(fp,"V2 x: %f,y: %f,z %f \n",TI.V2.x,TI.V2.y,TI.V2.z);
			*/

			E.v1 = PointInTriangle(TI.T2,TI.V1);
			E.v2 = PointInTriangle(TI.T2,TI.V2);


			if(!E.v1.valid()){
				context->OutputToConsoleEx("x: %f,y: %f,z %f",TI.V1.x,TI.V1.y,TI.V1.z);
				context->OutputToConsoleEx("x: %f,y: %f,z %f",E.v1.pointInTriangle.x,E.v1.pointInTriangle.y,E.v1.pointInTriangle.z);
				DEBUGBREAK
			}

			if(!E.v2.valid()){
				context->OutputToConsoleEx("x: %f,y: %f,z %f",TI.V2.x,TI.V2.y,TI.V2.z);
				context->OutputToConsoleEx("x: %f,y: %f,z %f",E.v2.pointInTriangle.x,E.v2.pointInTriangle.y,E.v2.pointInTriangle.z);
				DEBUGBREAK
			}


			//context->OutputToConsoleEx("T2: %d",TI.T2.faceIndex);
			__ShowVxVector(context,E.v1.pointInTriangle);
			__ShowVxVector(context,E.v1.point);
			__ShowVxVector(context,E.v2.pointInTriangle);
			__ShowVxVector(context,E.v2.point);


			T2Polygons[TI.T2.faceIndex].Edges.push_back(E);
		}
		
	}//for i

	//context->OutputToConsoleEx("T1Triangles.size %d",T1Triangles.size());
	//context->OutputToConsoleEx("T2Triangles.size %d",T2Triangles.size());

	std::vector<std::vector<Triangle> > incrTriangles;

	
	for(std::set<int>::iterator it = T1Triangles.begin();it != T1Triangles.end();it++){
		//context->OutputToConsoleEx("T1 %d ~~~~~~~~~~~~~~~~~~~~",*it);
		T1Polygons[*it].MakeChains(context);
		for(std::vector<Chain>::iterator ChainIter = T1Polygons[*it].chains.begin();
			ChainIter!=T1Polygons[*it].chains.end();ChainIter++){

				/*
				for(std::list<Edge>::iterator EdgeIter = ChainIter->Edges.begin();
					EdgeIter != ChainIter->Edges.end();EdgeIter++){
						context->OutputToConsoleEx("E.v1(%f,%f,%f) E.v2(%f,%f,%f)",EdgeIter->v1.pointInTriangle.v[0],
							EdgeIter->v1.pointInTriangle.v[1],EdgeIter->v1.pointInTriangle.v[2],EdgeIter->v2.pointInTriangle.v[0],
							EdgeIter->v2.pointInTriangle.v[1],EdgeIter->v2.pointInTriangle.v[2]);
				}
				*/
				
				

		}//for ChainIter
		std::vector<Area> Areas = T1Polygons[*it].SortChains(context);
		
		//context->OutputToConsoleEx("Areas.size() %d",Areas.size());
		//context->OutputToConsoleEx("Chains %d",T1Polygons[*it].chains.size());
		for(std::vector<Chain>::iterator ChainIter = T1Polygons[*it].chains.begin();
			ChainIter!=T1Polygons[*it].chains.end();ChainIter++){

				/*
				for(std::list<Edge>::iterator EdgeIter = ChainIter->Edges.begin();
					EdgeIter != ChainIter->Edges.end();EdgeIter++){
						context->OutputToConsoleEx("E.v1(%f,%f,%f) E.v2(%f,%f,%f)",EdgeIter->v1.pointInTriangle.v[0],
							EdgeIter->v1.pointInTriangle.v[1],EdgeIter->v1.pointInTriangle.v[2],EdgeIter->v2.pointInTriangle.v[0],
							EdgeIter->v2.pointInTriangle.v[1],EdgeIter->v2.pointInTriangle.v[2]);
				}
				*/
				

		}//for ChainIter
        
		for(size_t i = 0;i < Areas.size();i++){
            std::vector<Triangle> triangles = Areas[i].Triangulation(context);

		    //context->OutputToConsoleEx("Triangles %d",triangles.size());
			incrTriangles.push_back(triangles);
		}

	}
	

	for(std::set<int>::iterator it = T2Triangles.begin();it != T2Triangles.end();it++){
		//context->OutputToConsoleEx("T2 %d ~~~~~~~~~~~~~~~~~~~~",*it);

		//context->OutputToConsoleEx("Edges %d",T2Polygons[*it].Edges.size());

		/*
		for(std::list<Edge>::iterator EdgeIter = T2Polygons[*it].Edges.begin();
			EdgeIter != T2Polygons[*it].Edges.end();EdgeIter++){
			context->OutputToConsoleEx("E.v1(%f,%f,%f) E.v2(%f,%f,%f)",EdgeIter->v1.pointInTriangle.v[0],
			EdgeIter->v1.pointInTriangle.v[1],EdgeIter->v1.pointInTriangle.v[2],EdgeIter->v2.pointInTriangle.v[0],
			EdgeIter->v2.pointInTriangle.v[1],EdgeIter->v2.pointInTriangle.v[2]);
		}
		*/
		
		T2Polygons[*it].MakeChains(context);
		//context->OutputToConsoleEx("Chains %d",T2Polygons[*it].chains.size());
		
		
		std::vector<Area> Areas = T2Polygons[*it].SortChains(context);
		
		/*
		context->OutputToConsoleEx("Areas.size() %d",Areas.size());
		context->OutputToConsoleEx("Chains %d",T2Polygons[*it].chains.size());
		for(std::vector<Chain>::iterator ChainIter = T2Polygons[*it].chains.begin();
			ChainIter!=T2Polygons[*it].chains.end();ChainIter++){

				
				for(std::list<Edge>::iterator EdgeIter = ChainIter->Edges.begin();
					EdgeIter != ChainIter->Edges.end();EdgeIter++){
						context->OutputToConsoleEx("E.v1(%f,%f,%f) E.v2(%f,%f,%f)",EdgeIter->v1.pointInTriangle.v[0],
							EdgeIter->v1.pointInTriangle.v[1],EdgeIter->v1.pointInTriangle.v[2],EdgeIter->v2.pointInTriangle.v[0],
							EdgeIter->v2.pointInTriangle.v[1],EdgeIter->v2.pointInTriangle.v[2]);
				}
				
		}//for ChainIter
		*/
		
		
		for(size_t i = 0;i < Areas.size();i++){
            std::vector<Triangle> triangles = Areas[i].Triangulation(context);
			// this is mesh2 in mesh1
			//so we should flip triangles

			for(std::vector<Triangle>::iterator tri = triangles.begin();tri != triangles.end(); tri++){
				std::swap(tri->v[1],tri->v[2]);
			}

		    //context->OutputToConsoleEx("Triangles %d",triangles.size());
			incrTriangles.push_back(triangles);
		}
		
		
	}


	CKVINDEX* faceIndex1Back = new CKVINDEX[Mesh1FaceCount * 3];

	//memcpy(faceIndex1Back,faceIndex1,sizeof(CKVINDEX) * Mesh1FaceCount * 3);

	memset(faceIndex1Back,0,sizeof(CKVINDEX) * Mesh1FaceCount* 3);

	int remainingFace = 0;

	int i = 0;

	/*
	for(std::vector<bool>::iterator it = Mesh1TriangleMarks.begin();it != Mesh1TriangleMarks.end();it++,i++){
		if(!*it && (Mesh1InMesh2[faceIndex1[i * 3]] != In
			&& Mesh1InMesh2[faceIndex1[i * 3 + 1]] != In
			&& Mesh1InMesh2[faceIndex1[i * 3 + 2]] != In)){
			faceIndex1Back[remainingFace * 3] = faceIndex1[i * 3];
			faceIndex1Back[remainingFace * 3 + 1] = faceIndex1[i * 3 + 1];
			faceIndex1Back[remainingFace * 3 + 2] = faceIndex1[i * 3 + 2];
			remainingFace++;
		}
	}
	*/

	for(std::vector<bool>::iterator it = Mesh1TriangleMarks.begin();it != Mesh1TriangleMarks.end();it++,i++){
		if(!*it && (Mesh1InMesh2[faceIndex1[i * 3]] == Out
			&& Mesh1InMesh2[faceIndex1[i * 3 + 1]] == Out
			&& Mesh1InMesh2[faceIndex1[i * 3 + 2]] == Out)){
			faceIndex1Back[remainingFace * 3] = faceIndex1[i * 3];
			faceIndex1Back[remainingFace * 3 + 1] = faceIndex1[i * 3 + 1];
			faceIndex1Back[remainingFace * 3 + 2] = faceIndex1[i * 3 + 2];
			remainingFace++;
		}
	}

	

	mesh1->SetFaceCount(remainingFace);
	faceIndex1 = mesh1->GetFacesIndices();

	memcpy(faceIndex1,faceIndex1Back,remainingFace * 3 * sizeof(CKVINDEX));

	i = 0;
	for(std::vector<bool>::iterator it = Mesh2TriangleMarks.begin();it != Mesh2TriangleMarks.end();it++,i++){
		if(!*it && Mesh2InMesh1[faceIndex2[i * 3]] == In
			&& Mesh2InMesh1[faceIndex2[i * 3 + 1]] == In
			&& Mesh2InMesh1[faceIndex2[i * 3 + 2]] == In){
			NewPoints.push_back( (*((VxVector*)(Mesh2Vertices + faceIndex2[i * 3] * Mesh2vStride )) + displacement) / scale1);
			
			//swap!
			NewPoints.push_back( (*((VxVector*)(Mesh2Vertices + faceIndex2[i * 3 + 2] * Mesh2vStride )) + displacement) / scale1);
			NewPoints.push_back( (*((VxVector*)(Mesh2Vertices + faceIndex2[i * 3 + 1] * Mesh2vStride )) + displacement) / scale1);
		}
	}

	for(std::vector<std::vector<Triangle> >::iterator it = incrTriangles.begin(); it != incrTriangles.end();it++){
		for(std::vector<Triangle>::iterator itit = it->begin();itit != it->end();itit++){
			itit->v[0] /= scale1;
			itit->v[1] /= scale1;
			itit->v[2] /= scale1;

			NewPoints.push_back(itit->v[0]);
			NewPoints.push_back(itit->v[1]);
			NewPoints.push_back(itit->v[2]);
			/*
			
			context->OutputToConsoleEx("~~~~~~~~~~~~~");

			ShowVxVector(context,itit->v[0]);
			ShowVxVector(context,itit->v[1]);
			ShowVxVector(context,itit->v[2]);

			*/
			
		}
	}

	int NewPointsCount = Mesh1VCount + NewPoints.size();

	mesh1->SetVertexCount(NewPointsCount);


	for(int i = Mesh1VCount;i < NewPointsCount;i++){
		mesh1->SetVertexPosition(i,&NewPoints[i - Mesh1VCount]);
	}

	mesh1->ModifierVertexMove(TRUE,TRUE);

	int NewIndexsCount = remainingFace + NewPoints.size()/3;

	mesh1->SetFaceCount(NewIndexsCount);

	faceIndex1 = mesh1->GetFacesIndices();

	for(int i = remainingFace;i < NewIndexsCount;i++){
		faceIndex1[i * 3] = Mesh1VCount + (i - remainingFace) * 3;
		faceIndex1[i * 3 + 1] = Mesh1VCount + (i - remainingFace) * 3 + 1;
		faceIndex1[i * 3 + 2] = Mesh1VCount + (i - remainingFace) * 3 + 2;
	}

	for(int i = 0;i < NewIndexsCount;i++){
		mesh1->SetFaceMaterial(i,mesh1->GetMaterial(0));
	}


	delete[] faceIndex1Back;
	
	mesh1->UVChanged();

	//context->OutputToConsoleEx("VerTex Count %d",mesh1->GetVertexCount());
	mesh1->BuildNormals();
	mesh1->BuildFaceNormals();
	mesh1->NormalChanged();


	/*
	for(int i = 0;i < NewPointsCount;i++){
		VxVector norm;
		mesh1->GetVertexNormal(i,&norm);
		ShowVxVector(context,norm);
		ShowVxVector(context,mesh1->GetFaceNormal(i));
		context->OutputToConsoleEx("i: %d , color %d",i,mesh1->GetVertexSpecularColor(i));
	}
	*/


}


int GetAllVertexBuildingBlock(const CKBehaviorContext& BehContext)
{

	

	CKBehavior* beh = BehContext.Behavior;

	CKContext* context = beh->GetCKContext();
	//context->OutputToConsoleEx("before all");

    CK3dEntity *ent = (CK3dEntity *)beh->GetTarget();

	VxVector myScale;

	//myScale = VxVector(1.0f,1.0f,1.0f);

	ent->GetScale(&myScale);

	CK3dEntity* watchEnt = (CK3dEntity*)beh->GetInputParameterObject(0);

	VxBbox box1 = ent->GetBoundingBox();
	VxBbox box2 = watchEnt->GetBoundingBox();

	VxVector watchScale;

	watchEnt->GetScale(&watchScale);

	//watchScale = VxVector(1.0f,1.0f,1.0f);

	VxVector watchPosition;
	watchEnt->GetPosition(&watchPosition);

	context->OutputToConsoleEx("watchPosition x: %f,y: %f,z %f",watchPosition.x,watchPosition.y,watchPosition.z);

	VxVector entPosition;
	ent->GetPosition(&entPosition);

	context->OutputToConsoleEx("entPosition x: %f,y: %f,z %f",entPosition.x,entPosition.y,entPosition.z);

	VxVector displacement = watchPosition - entPosition;

	context->OutputToConsoleEx("displacement x: %f,y: %f,z %f",displacement.x,displacement.y,displacement.z);

	if(!ent) return CKBR_OWNERERROR;

	CKMesh *mesh = ent->GetCurrentMesh();

	CKMesh *WatchMesh = (CKMesh *)watchEnt->GetMesh(0);

	if(!mesh){
		context->OutputToConsole("can't fetch mesh",TRUE);
		return CKBR_OK;
	}


	if(!WatchMesh){
		context->OutputToConsole("can't fetch WatchMesh",TRUE);
		return CKBR_OK;
	}

	//std::pair<std::vector<int>,std::vector<bool> > MeshInWatchMesh = Mesh1PointsInMesh2(mesh,WatchMesh,-watchPosition,context);
	//std::pair<std::vector<int>,std::vector<bool> > WatchMeshInMesh = Mesh1PointsInMesh2(WatchMesh,mesh, watchPosition,context);

	try{
	    //std::vector<TriangleIntersection> TriangleIntersection = IntersectedTrianglesOf2Mesh(mesh,WatchMesh,watchPosition,context);
	    CutMesh1ByMesh2(mesh,WatchMesh,box1,box2,displacement,myScale,watchScale,context);
	}catch(std::string str){
		throw str;
	}
	
	beh->ActivateInput(0,FALSE);
	return CKBR_OK;
}

int GetAllVertexBuildingBlockCallBack(const CKBehaviorContext& BehContext)
{
	switch (BehContext.CallbackMessage)
	{
		case CKM_BEHAVIORATTACH:
			break;
		case CKM_BEHAVIORDETACH:
			break;
		case CKM_BEHAVIORDELETE:
			break;
		case CKM_BEHAVIOREDITED:
			break;
		case CKM_BEHAVIORSETTINGSEDITED:
			break;
		case CKM_BEHAVIORLOAD:
			break;
		case CKM_BEHAVIORPRESAVE:
			break;
		case CKM_BEHAVIORPOSTSAVE:
			break;
		case CKM_BEHAVIORRESUME:
			break;
		case CKM_BEHAVIORPAUSE:
			break;
		case CKM_BEHAVIORRESET:
			break;
		case CKM_BEHAVIORNEWSCENE:
			break;
		case CKM_BEHAVIORDEACTIVATESCRIPT:
			break;
		case CKM_BEHAVIORACTIVATESCRIPT:
			break;
		case CKM_BEHAVIORREADSTATE:
			break;

	}
	return CKBR_OK;
}


