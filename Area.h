#ifndef __AREA_H__
#define __AREA_H__

#include "Chain.h"

//Fix Me: release pointers

//use RAII
template <typename T>
class Memorydeleter{
public:
	Memorydeleter(){}
	~Memorydeleter(){
		for(size_t i = 0;i<vec.size();i++){
			delete vec[i];
		}
	}
	void add(T ptr){
		vec.push_back(ptr);
	}
private:
	std::vector<T> vec;
};

class Area{
public:
	Chain Boundary;
	std::vector<Chain> innerBoundarys;
	VxVector ParentNorm;
	std::vector<Triangle> Triangulation(CKContext* context){

		std::vector<Triangle> vec;

		ChainNode* root = 0;
		ChainNode* curr = 0;
		ChainNode* prev = 0;

		Memorydeleter<ChainNode*> deleter;

		int remainingEdges = Boundary.Edges.size();
		//just return a triangle
		
		//context->OutputToConsoleEx("begin Area");

		for(std::list<Edge>::iterator it = Boundary.Edges.begin(); it != Boundary.Edges.end();it++){
			curr = new ChainNode();
			deleter.add(curr);
			curr->V = (*it).v1;

			//context->OutputToConsoleEx("begin a node");

			ShowVxVector(context,(*it).v1.pointInTriangle);
			ShowVxVector(context,(*it).v2.pointInTriangle);

			//context->OutputToConsoleEx("root %p, curr %p, prev %p",root,curr,prev);
			if(!root){
				root = curr;
			}else{
				prev->next = curr;
				curr->prev = prev;
			}
			prev = curr;
		}

		curr->next = root;
		root->prev = curr;

		ChainNode* checkptr = root;
		do{
			VxVector& p1 = checkptr->prev->V.point;
			VxVector& p2 = checkptr->V.point;
			VxVector& p3 = checkptr->next->V.point;

			VxVector l1 = p2 - p1;
			VxVector l2 = p3 - p2;

			VxVector xtemp = VxVectorCrossProduct(l1,l2);
				

			//if(AlmostEqualZero(xtemp.x) && AlmostEqualZero(xtemp.y) && AlmostEqualZero(xtemp.z)){
			if(xtemp.Magnitude() < 0.00001){

				//on one line;
				//remove Checkptr
				//context->OutputToConsoleEx("remove checkptr");
				checkptr->prev->next = checkptr->next;
				checkptr->next->prev = checkptr->prev;

				if(checkptr == root)
					root = checkptr->prev;
				remainingEdges--;
				//the triangle is too small
				if(remainingEdges < 3)
			           return vec;
			}

		    checkptr = checkptr->next;
		}while(checkptr!=root);

		
		bool GotOne = false;
		
		ChainNode* checkFatalError = root;
		do{
			VxVector& p1 = checkFatalError->prev->V.point;
			VxVector& p2 = checkFatalError->V.point;
			VxVector& p3 = checkFatalError->next->V.point;

			VxVector p2norm = VxVectorCrossProduct( p1 - p2,p3 - p2);
			p2norm.Normalize();

			if(SameVertex(p2norm,ParentNorm)){
				//context->OutputToConsoleEx("got one");
				//context->OutputToConsoleEx("p2norm (%f,%f,%f) ParentNorm (%f,%f,%f)",p2norm.x,p2norm.y,p2norm.z,ParentNorm.x,ParentNorm.y,ParentNorm.z);
				GotOne = true;
				break;
			}


			checkFatalError = checkFatalError->next;
		}while(checkFatalError != root);

		if(!GotOne){
			context->OutputToConsoleEx("flip the chain");
			checkFatalError = root;
			do{
				ChainNode* nextNode = checkFatalError->next;
				std::swap(checkFatalError->next,checkFatalError->prev);
				checkFatalError = nextNode;
			}while(checkFatalError != root);
		}
		
		

		int reject = 0;

		if(remainingEdges < 3)
			DEBUGBREAK

		while(true){
			if(remainingEdges == 3){
			
				

				Triangle ret;
				ret.v[0] = root->V.point;
				ret.v[1] = root->next->V.point;
				ret.v[2] = root->next->next->V.point;

				VxVector& p1 = root->prev->V.point;
				VxVector& p2 = root->V.point;
				VxVector& p3 = root->next->V.point;

				VxVector p2norm = VxVectorCrossProduct( p1 - p2,p3 - p2);
				p2norm.Normalize();

				if(!SameVertex(p2norm,ParentNorm)){
					std::swap(ret.v[1],ret.v[2]);
				}


				vec.push_back(ret);
				return vec;
			}
			else{
				
				ChainNode* p = root;
				do{
					VxVector& p1 = p->prev->V.point;
					VxVector& p2 = p->V.point;
					VxVector& p3 = p->next->V.point;

					VxVector p2norm = VxVectorCrossProduct( p1 - p2,p3 - p2);
					VxVector t = p2norm;
					p2norm.Normalize();

					if(reject > remainingEdges + 1){
						context->OutputToConsoleEx("Fatal Error");
						return vec;
					}

					if(!SameVertex(p2norm,ParentNorm)){
						//context->OutputToConsoleEx("reject by norm");
						//context->OutputToConsoleEx("t (%f,%f,%f)",t.x,t.y,t.z);
						//context->OutputToConsoleEx("p2norm (%f,%f,%f) ParentNorm (%f,%f,%f)",p2norm.x,p2norm.y,p2norm.z,ParentNorm.x,ParentNorm.y,ParentNorm.z);
						reject++;
						p = p->next;
						continue;
					}

					Triangle test(p1,p2,p3,0);

					bool valid = true;

					//check if valid
					ChainNode* inner = root;
					do{
						
						if(inner == p || inner == p->prev || inner == p->next){
							inner = inner->next;
							continue;
						}
						VxVector pointInTriangle = test.PositionOfPointInTriangle(inner->V.point);
						if(InTheRange(0,1,pointInTriangle.x) && InTheRange(0,1,pointInTriangle.y) && InTheRange(0,1,pointInTriangle.z)){
							//point in triangle
							valid = false;
							break;
						}
						inner = inner->next;
					}while(inner != root);

					if(valid){
                        Triangle temp;
				        temp.v[0] = p1;
				        temp.v[1] = p2;
				        temp.v[2] = p3;
					    vec.push_back(temp);
						p->prev->next = p->next;
						p->next->prev = p->prev;
						root = p->next;
						remainingEdges--;
						reject = 0;
						//return vec;
						break;
					}

					

					p = p->next;
				}while(p!=root);
			}//else
		}

		return vec;

	}
};

#endif