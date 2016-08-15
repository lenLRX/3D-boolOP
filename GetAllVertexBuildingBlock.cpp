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
	CKObjectDeclaration *od = CreateCKObjectDeclaration("GetAllVertexBuildingBlock");
	
	od->SetType(CKDLL_BEHAVIORPROTOTYPE);
	od->SetVersion(0x00000001);
	od->SetCreationFunction(CreateGetAllVertexBuildingBlockProto);
	od->SetDescription("Enter your description here");
	od->SetCategory("UserBBs");
	od->SetGuid(CKGUID(0x4acca9a7, 0xb22da652));
	od->SetAuthorGuid(CKGUID(0x56495254,0x4f4f4c53));
	od->SetAuthorName("Virtools");
	od->SetCompatibleClassId(CKCID_BEOBJECT);

	return od;
}

CKERROR CreateGetAllVertexBuildingBlockProto(CKBehaviorPrototype** pproto)
{
	CKBehaviorPrototype *proto = CreateCKBehaviorPrototype("GetAllVertexBuildingBlock");
	if (!proto) {
		return CKERR_OUTOFMEMORY;
	}

//---     Inputs declaration
	proto->DeclareInput("In0");

//---     Outputs declaration
	proto->DeclareOutput("Out0");


	proto->DeclareInParameter("WatchObject", CKPGUID_3DENTITY );

//----	Local Parameters Declaration

//----	Settings Declaration
	proto->SetBehaviorFlags((CK_BEHAVIOR_FLAGS)(CKBEHAVIOR_TARGETABLE));
	proto->SetBehaviorCallbackFct(GetAllVertexBuildingBlockCallBack, CKCB_BEHAVIORATTACH|CKCB_BEHAVIORDETACH|CKCB_BEHAVIORDELETE|CKCB_BEHAVIOREDITED|CKCB_BEHAVIORSETTINGSEDITED|CKCB_BEHAVIORLOAD|CKCB_BEHAVIORPRESAVE|CKCB_BEHAVIORPOSTSAVE|CKCB_BEHAVIORRESUME|CKCB_BEHAVIORPAUSE|CKCB_BEHAVIORRESET|CKCB_BEHAVIORRESET|CKCB_BEHAVIORDEACTIVATESCRIPT|CKCB_BEHAVIORACTIVATESCRIPT|CKCB_BEHAVIORREADSTATE, NULL);
	proto->SetFunction(GetAllVertexBuildingBlock);

	*pproto = proto;
	return CK_OK;
}

int PointInBody(VxVector pt,VxVector scale,CKMesh* mesh,int flag,CKContext* context){

	CKDWORD vStride=0;

    const BYTE* vxv = (BYTE*)mesh->GetModifierVertices(&vStride);
	int vcount = mesh->GetModifierVertexCount();

	int faceCount = mesh->GetFaceCount();
	
	const CKVINDEX* faceIndex = mesh->GetFacesIndices();

	VxVector line = VxVector(-1,-1,-1);

	line = VxVector(-(float)(rand() + 1)/(float)RAND_MAX,-(float)(rand() + 1)/(float)RAND_MAX,-(float)(rand() + 1)/(float)RAND_MAX);

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

		
		
		

		bool b = triangle.VxVevtorIntersectTriangle(pt,line,flag,intersectionPoint);

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
		

		if(b){
			IntersectionCount++;
	    }
	}

	
	return IntersectionCount;
}

std::pair<std::vector<int>,std::vector<bool> > Mesh1PointsInMesh2(CKMesh* mesh1,CKMesh* mesh2,VxVector displacement,VxVector scale1,VxVector scale2,int flag,CKContext* context){
	CKDWORD Mesh1vStride=0;
	CKDWORD Mesh2vStride=0;
	const BYTE* Mesh1Vertices = (BYTE*)mesh1->GetModifierVertices(&Mesh1vStride);
	const BYTE* Mesh2Vertices = (BYTE*)mesh2->GetModifierVertices(&Mesh2vStride);
	int Mesh1VCount = mesh1->GetModifierVertexCount();
	int Mesh2VCount = mesh2->GetModifierVertexCount();

	std::pair<std::vector<int>,std::vector<bool> > ret;
	int count = 0;
	
	
	for(int i=0; i<Mesh1VCount ;i++,Mesh1Vertices+=Mesh1vStride) {

		VxVector TransformedPos = *(VxVector*)Mesh1Vertices;
		for(int j = 0;j < 3;j++){
			TransformedPos.v[j] *= scale1.v[j];
		}

		TransformedPos += displacement;

		int IntersectionCount = PointInBody(TransformedPos,scale2,mesh2,flag,context);

		if(IntersectionCount % 2 == 1){
			ret.first.push_back(i);
			ret.second.push_back(true);
		}
		else
			ret.second.push_back(false);
	}
	

	return ret;
}

std::vector<TriangleIntersection> IntersectedTrianglesOf2Mesh(CKMesh* mesh1,CKMesh* mesh2,
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

		for(int k = 0;k < 3;k++){
			T1v1.v[k] *= scale1.v[k];
			T1v2.v[k] *= scale1.v[k];
			T1v3.v[k] *= scale1.v[k];
		}


		Triangle T1(T1v1,T1v2,T1v3,i);

		for(int j = 0;j < Mesh2FaceCount;j++){
	        CKVINDEX T2f1 = *(faceIndex2 + 3 * j);
			CKVINDEX T2f2 = *(faceIndex2 + 3 * j + 1);
			CKVINDEX T2f3 = *(faceIndex2 + 3 * j + 2);

			//context->OutputToConsoleEx("mesh2 faces %d,%d,%d",T2f1,T2f2,T2f3);

			VxVector T2v1 = *((VxVector*)(Mesh2Vertices + T2f1*Mesh2vStride )) + displacement;
			VxVector T2v2 = *((VxVector*)(Mesh2Vertices + T2f2*Mesh2vStride )) + displacement;
			VxVector T2v3 = *((VxVector*)(Mesh2Vertices + T2f3*Mesh2vStride )) + displacement;

			for(int k = 0;k < 3;k++){
			    T2v1.v[k] *= scale2.v[k];
			    T2v2.v[k] *= scale2.v[k];
			    T2v3.v[k] *= scale2.v[k];
		    }

			Triangle T2(T2v1,T2v2,T2v3,j);

			//TODO:check if it is in triangle

			Intersections intersections = IntersectWith(context,T1,T2);
			if(intersections.size()){
				IntersectionCount += intersections.size();
				
				for(size_t idx = 0;idx < intersections.size();idx++){
					/*
					PointInTriangle p1(T1,intersections[idx].V1);
					if(!p1.valid()){
						__ShowVxVector(context,p1.pointInTriangle);
						DEBUGBREAK
					}

					PointInTriangle p2(T1,intersections[idx].V2);
					if(!p2.valid()){
						__ShowVxVector(context,p2.pointInTriangle);
						DEBUGBREAK
					}

					PointInTriangle p3(T2,intersections[idx].V1);
					if(!p3.valid()){
						__ShowVxVector(context,p3.pointInTriangle);
						DEBUGBREAK
					}

					PointInTriangle p4(T2,intersections[idx].V1);
					if(!p4.valid()){
						__ShowVxVector(context,p4.pointInTriangle);
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

void CutMesh1ByMesh2(CKMesh* mesh1,CKMesh* mesh2,
					 VxVector displacement,VxVector scale1,VxVector scale2,CKContext* context){
    srand(time(NULL));
    std::pair<std::vector<int>,std::vector<bool> > Mesh1InMesh2 = Mesh1PointsInMesh2(mesh1,mesh2,-displacement,scale1,scale2,1,context);
    std::pair<std::vector<int>,std::vector<bool> > Mesh2InMesh1 = Mesh1PointsInMesh2(mesh2,mesh1, displacement,scale2,scale1,-1,context);

	size_t s = Mesh2InMesh1.second.size();
	for(size_t i = 0;i < s;i++){
		// visibility is inverse for Mesh2InMesh1
		Mesh2InMesh1.second[i] = !Mesh2InMesh1.second[i];
	}

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

	std::vector<TriangleIntersection> TriangleIntersections = IntersectedTrianglesOf2Mesh(mesh1,mesh2,displacement,scale1,scale2,context);

	size_t IntersectionSize = TriangleIntersections.size();

	std::vector<VxVector> NewVertices;
	int NewVerticesIndex = 0;
	std::map<VxVector,int,VxVectorLess> VxVector2Index;

	for(int i = 0 ; i < Mesh1VCount; i++){
		if(Mesh1InMesh2.second[i] == false){
			VxVector temp(*((VxVector*)(Mesh1Vertices + i * Mesh1vStride)));
			if(VxVector2Index.find(temp) == VxVector2Index.end()){
				NewVertices.push_back(temp);
			    VxVector2Index[temp] = NewVerticesIndex;
			    NewVerticesIndex++;
			}
		}
	}

	for(int i = 0 ; i < Mesh2VCount; i++){
		if(Mesh2InMesh1.second[i] == true){
			VxVector temp(*((VxVector*)(Mesh2Vertices + i * Mesh2vStride)));
			if(VxVector2Index.find(temp) == VxVector2Index.end()){
				NewVertices.push_back(temp);
			    VxVector2Index[temp] = NewVerticesIndex;
			    NewVerticesIndex++;
			}
		}
	}

	std::vector<bool> Mesh1TriangleMarks(Mesh1FaceCount,false);
	std::vector<bool> Mesh2TriangleMarks(Mesh2FaceCount,false);

	std::vector<int> NewFaces;

	std::vector<Polygon> T1Polygons(Mesh1FaceCount);
	std::set<int> T1Triangles;
	std::vector<Polygon> T2Polygons(Mesh2FaceCount);
	std::set<int> T2Triangles;

	for(size_t i = 0; i < IntersectionSize ; i++){
		TriangleIntersection& TI = TriangleIntersections[i];
		Mesh1TriangleMarks[TI.T1.faceIndex] = true;
		Mesh2TriangleMarks[TI.T2.faceIndex] = true;
		if(VxVector2Index.find(TI.V1) == VxVector2Index.end()){
		    NewVertices.push_back(TI.V1);
			VxVector2Index[TI.V1] = NewVerticesIndex;
			NewVerticesIndex++;
		}

		if(VxVector2Index.find(TI.V2) == VxVector2Index.end()){
		    NewVertices.push_back(TI.V2);
			VxVector2Index[TI.V2] = NewVerticesIndex;
			NewVerticesIndex++;
		}

		Edge E;

		if(TI.T1valid){
			T1Polygons[TI.T1.faceIndex].triangle = TI.T1;

			T1Triangles.insert(TI.T1.faceIndex);

			for(int j = 0; j < 3;j++){
				T1Polygons[TI.T1.faceIndex].triangle.visible[j]
				= !Mesh1InMesh2.second[*(faceIndex1 + 3 * TI.T1.faceIndex + j)];
			}

			E.v1 = PointInTriangle(TI.T1,TI.V1);
			E.v2 = PointInTriangle(TI.T1,TI.V2);

			T1Polygons[TI.T1.faceIndex].Edges.push_back(E);
		}

		if(TI.T2valid){

			T2Polygons[TI.T2.faceIndex].triangle = TI.T2;

			T2Triangles.insert(TI.T2.faceIndex);

			for(int j = 0; j < 3;j++){
				T2Polygons[TI.T2.faceIndex].triangle.visible[j] 
				= !Mesh2InMesh1.second[*(faceIndex2 + 3 * TI.T2.faceIndex + j)];
			}
			
			E.v1 = PointInTriangle(TI.T2,TI.V1);
			E.v2 = PointInTriangle(TI.T2,TI.V2);

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

	for(std::vector<bool>::iterator it = Mesh1TriangleMarks.begin();it != Mesh1TriangleMarks.end();it++,i++){
		if(!*it && Mesh1InMesh2.second[faceIndex1[i * 3]] == false
			&& Mesh1InMesh2.second[faceIndex1[i * 3 + 1]] == false
			&& Mesh1InMesh2.second[faceIndex1[i * 3 + 2]] == false){
			faceIndex1Back[remainingFace * 3] = faceIndex1[i * 3];
			faceIndex1Back[remainingFace * 3 + 1] = faceIndex1[i * 3 + 1];
			faceIndex1Back[remainingFace * 3 + 2] = faceIndex1[i * 3 + 2];
			remainingFace++;
		}
	}

	

	mesh1->SetFaceCount(remainingFace);
	faceIndex1 = mesh1->GetFacesIndices();

	memcpy(faceIndex1,faceIndex1Back,remainingFace * 3 * sizeof(CKVINDEX));

	std::vector<VxVector> NewPoints;

	i = 0;
	for(std::vector<bool>::iterator it = Mesh2TriangleMarks.begin();it != Mesh2TriangleMarks.end();it++,i++){
		//notice: we had fliped it before
		if(!*it && Mesh2InMesh1.second[faceIndex2[i * 3]] == false
			&& Mesh2InMesh1.second[faceIndex2[i * 3 + 1]] == false
			&& Mesh2InMesh1.second[faceIndex2[i * 3 + 2]] == false){
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
	    CutMesh1ByMesh2(mesh,WatchMesh,displacement,myScale,watchScale,context);
	}catch(std::string str){
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


