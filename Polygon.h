#ifndef __POLYGON_H__
#define __POLYGON_H__

#include <vector>
#include <list>
#include <algorithm>
#include "MyUtility.h"
#include "Edge.h"
#include "Triangle.h"
#include "Chain.h"
#include "Area.h"

struct SortByGreater
{
	SortByGreater(int idx):i(idx){}
	int i;
	bool operator () (const Chain& lhs,const Chain& rhs)const {
		return lhs.max[i] > rhs.max[i];
	}
};

static int additionMsg2Vertex[4]={-1,2,1,0};

class Polygon
{
public:
	Polygon(){}
	std::list<Edge> Edges;
	std::vector<Chain> chains;
	Triangle triangle;

	void MakeChains(CKContext* context){
		for(std::list<Edge>::iterator it1 = Edges.begin();it1!=Edges.end();it1++){
			for(std::list<Edge>::iterator it2 = Edges.begin();it2!=Edges.end();it2++){
				if(it1 == it2){
					continue;
				}

				if((SameVertex(it1->v1.point,it2->v1.point) && SameVertex(it1->v2.point,it2->v2.point))
					&& (SameVertex(it1->v1.point,it2->v2.point) && SameVertex(it1->v2.point,it2->v1.point))){
					it1->valid = false;
					it2->valid = false;
				}
			}
		}

		for(std::list<Edge>::iterator it = Edges.begin();it!=Edges.end();){
			if(!it->valid || SameVertex(it->v1.point,it->v2.point)){
				it = Edges.erase(it); 
			}else{
				it++;
			}
		}



		while(Edges.size()){
			Edge current = Edges.front();
			Chain _chain;
			_chain.Edges.push_back(current);


			for(int i = 0;i < 3;i++){
							if(current.v1.pointInTriangle.v[i] > _chain.max[i])
								_chain.max[i] = current.v1.pointInTriangle.v[i];
							if(current.v1.pointInTriangle.v[i] < _chain.min[i])
								_chain.min[i] = current.v1.pointInTriangle.v[i];

							if(current.v2.pointInTriangle.v[i] > _chain.max[i])
								_chain.max[i] = current.v1.pointInTriangle.v[i];
							if(current.v2.pointInTriangle.v[i] < _chain.min[i])
								_chain.min[i] = current.v1.pointInTriangle.v[i];
						}

			Edges.erase(Edges.begin());
			bool onceGotEdge;
			do
			{
				onceGotEdge = false;
				for(std::list<Edge>::iterator it = Edges.begin(); it != Edges.end();){
					if(SameVertex(_chain.Edges.back().v2.point,it->v1.point)){
						if(_chain.Edges.back().v2.OnTheEdge >= 0){// && it->v1.OnTheEdge >=0){
							//at Edge of Triangle, chain must be Splited, now we just skip it,it will be processed Later
							it++;
							continue;
						}
						current = *it;

						for(int i = 0;i < 3;i++){
							if(current.v1.pointInTriangle.v[i] > _chain.max[i])
								_chain.max[i] = current.v1.pointInTriangle.v[i];
							if(current.v1.pointInTriangle.v[i] < _chain.min[i])
								_chain.min[i] = current.v1.pointInTriangle.v[i];

							if(current.v2.pointInTriangle.v[i] > _chain.max[i])
								_chain.max[i] = current.v1.pointInTriangle.v[i];
							if(current.v2.pointInTriangle.v[i] < _chain.min[i])
								_chain.min[i] = current.v1.pointInTriangle.v[i];
						}

						//context->OutputToConsoleEx("_chain.Edges.back().v2.point,it->v1.point");
						ShowVxVector(context,_chain.Edges.back().v1.pointInTriangle);
						ShowVxVector(context,_chain.Edges.back().v2.pointInTriangle);
						ShowVxVector(context,it->v1.pointInTriangle);
						ShowVxVector(context,it->v2.pointInTriangle);

						_chain.Edges.push_back(current);

						it = Edges.erase(it);
						onceGotEdge = true;
					}
					else if(SameVertex(_chain.Edges.back().v2.point,it->v2.point)){
						if(_chain.Edges.back().v2.OnTheEdge >= 0){// && it->v2.OnTheEdge >=0){
							//at Edge of Triangle, chain must be Splited
							it++;
							continue;
						}
						current = *it;

						for(int i = 0;i < 3;i++){
							if(current.v1.pointInTriangle.v[i] > _chain.max[i])
								_chain.max[i] = current.v1.pointInTriangle.v[i];
							if(current.v1.pointInTriangle.v[i] < _chain.min[i])
								_chain.min[i] = current.v1.pointInTriangle.v[i];

							if(current.v2.pointInTriangle.v[i] > _chain.max[i])
								_chain.max[i] = current.v1.pointInTriangle.v[i];
							if(current.v2.pointInTriangle.v[i] < _chain.min[i])
								_chain.min[i] = current.v1.pointInTriangle.v[i];
						}

						std::swap(current.v1,current.v2);


						//context->OutputToConsoleEx("_chain.Edges.back().v2.point,it->v2.point");
						ShowVxVector(context,_chain.Edges.back().v1.pointInTriangle);
						ShowVxVector(context,_chain.Edges.back().v2.pointInTriangle);
						ShowVxVector(context,current.v1.pointInTriangle);
						ShowVxVector(context,current.v2.pointInTriangle);

						_chain.Edges.push_back(current);
						it = Edges.erase(it);
						onceGotEdge = true;
					}
					else if(SameVertex(_chain.Edges.front().v1.point,it->v1.point)){
						if(_chain.Edges.front().v1.OnTheEdge >= 0){// && it->v1.OnTheEdge >=0){
							//at Edge of Triangle, chain must be Splited
							it++;
							continue;
						}
						current = *it;

						for(int i = 0;i < 3;i++){
							if(current.v1.pointInTriangle.v[i] > _chain.max[i])
								_chain.max[i] = current.v1.pointInTriangle.v[i];
							if(current.v1.pointInTriangle.v[i] < _chain.min[i])
								_chain.min[i] = current.v1.pointInTriangle.v[i];

							if(current.v2.pointInTriangle.v[i] > _chain.max[i])
								_chain.max[i] = current.v1.pointInTriangle.v[i];
							if(current.v2.pointInTriangle.v[i] < _chain.min[i])
								_chain.min[i] = current.v1.pointInTriangle.v[i];
						}

						std::swap(current.v1,current.v2);


						//context->OutputToConsoleEx("_chain.Edges.front().v1.point,it->v1.point");
						ShowVxVector(context,_chain.Edges.front().v1.pointInTriangle);
						ShowVxVector(context,_chain.Edges.front().v2.pointInTriangle);
						ShowVxVector(context,current.v1.pointInTriangle);
						ShowVxVector(context,current.v2.pointInTriangle);

						_chain.Edges.push_front(current);
						it = Edges.erase(it);
						onceGotEdge = true;
					}
					else if(SameVertex(_chain.Edges.front().v1.point,it->v2.point)){
						if(_chain.Edges.front().v1.OnTheEdge >= 0){// && it->v2.OnTheEdge >=0){
							//at Edge of Triangle, chain must be Splited
							it++;
							continue;
						}
						current = *it;

						for(int i = 0;i < 3;i++){
							if(current.v1.pointInTriangle.v[i] > _chain.max[i])
								_chain.max[i] = current.v1.pointInTriangle.v[i];
							if(current.v1.pointInTriangle.v[i] < _chain.min[i])
								_chain.min[i] = current.v1.pointInTriangle.v[i];

							if(current.v2.pointInTriangle.v[i] > _chain.max[i])
								_chain.max[i] = current.v1.pointInTriangle.v[i];
							if(current.v2.pointInTriangle.v[i] < _chain.min[i])
								_chain.min[i] = current.v1.pointInTriangle.v[i];
						}

						//context->OutputToConsoleEx("_chain.Edges.front().v1.point,it->v2.point");
						ShowVxVector(context,_chain.Edges.front().v1.pointInTriangle);
						ShowVxVector(context,_chain.Edges.front().v2.pointInTriangle);
						ShowVxVector(context,current.v1.pointInTriangle);
						ShowVxVector(context,current.v2.pointInTriangle);

						_chain.Edges.push_front(current);


						it = Edges.erase(it);
						onceGotEdge = true;
					}
					else
						it++;
				}//for
			}while(onceGotEdge);
			chains.push_back(_chain);
		}//while Edges.size()

		for(std::vector<Chain>::iterator it = chains.begin();
			it != chains.end();it++){
			if(it->Edges.front().v1.OnTheEdge >= 0 && it->Edges.back().v2.OnTheEdge >= 0){
				//both end on the Edge
				if(it->Edges.front().v1.OnTheEdge == it->Edges.back().v2.OnTheEdge){
					//They are on the same Edge
					it->type = ChainsOnSameEdge;
					it->additionalMsg = it->Edges.front().v1.OnTheEdge;
					//context->OutputToConsoleEx("They are on the same Edge additionMsg %d",it->additionalMsg);
				}
				else{
					it->type = ChainsOppositeToVertex;
					it->additionalMsg = it->Edges.front().v1.OnTheEdge + it->Edges.back().v2.OnTheEdge;
					//context->OutputToConsoleEx("They are not on the same Edge additionMsg %d",it->additionalMsg);
				}
			}
			else if(it->Edges.front().v1.OnTheEdge < 0 && it->Edges.back().v2.OnTheEdge < 0){
				it->type = EnclosedChains;
				//context->OutputToConsoleEx("They are EnClosedChains");
			}
			else{
				
				context->OutputToConsoleEx("size: %d", (int)it->Edges.size());
				context->OutputToConsoleEx("origin E.v1(%f,%f,%f) E.v2(%f,%f,%f)",it->Edges.front().v1.point.v[0],
							it->Edges.front().v1.point.v[1],it->Edges.front().v1.point.v[2],it->Edges.back().v2.point.v[0],
							it->Edges.back().v2.point.v[1],it->Edges.back().v2.point.v[2]);
				context->OutputToConsoleEx("E.v1(%f,%f,%f) E.v2(%f,%f,%f)",it->Edges.front().v1.pointInTriangle.v[0],
							it->Edges.front().v1.pointInTriangle.v[1],it->Edges.front().v1.pointInTriangle.v[2],it->Edges.back().v2.pointInTriangle.v[0],
							it->Edges.back().v2.pointInTriangle.v[1],it->Edges.back().v2.pointInTriangle.v[2]);
				
				//throw std::exception("OnTheEdge must be both >= 0 or < 0 !!");
				float min1 = 1.5f;
				float min2 = 1.5f;
				for(int i = 0;i < 3;i++){
					if(it->Edges.front().v1.pointInTriangle.v[i] < min1){
						min1 = it->Edges.front().v1.pointInTriangle.v[i];
						it->Edges.front().v1.OnTheEdge = i;
					}

					if(it->Edges.back().v2.pointInTriangle.v[i] < min2){
						min2 = it->Edges.back().v2.pointInTriangle.v[i];
						it->Edges.back().v2.OnTheEdge = i;
					}
				}

				if(it->Edges.front().v1.OnTheEdge >= 0 && it->Edges.back().v2.OnTheEdge >= 0){
					//both end on the Edge
					if(it->Edges.front().v1.OnTheEdge == it->Edges.back().v2.OnTheEdge){
						//They are on the same Edge
						it->type = ChainsOnSameEdge;
						it->additionalMsg = it->Edges.front().v1.OnTheEdge;
						//context->OutputToConsoleEx("They are on the same Edge additionMsg %d",it->additionalMsg);
					}
					else{
						it->type = ChainsOppositeToVertex;
						it->additionalMsg = it->Edges.front().v1.OnTheEdge + it->Edges.back().v2.OnTheEdge;
						//context->OutputToConsoleEx("They are not on the same Edge additionMsg %d",it->additionalMsg);
					}
				}
			}
			
		}
	}


	void check(){
		for(std::vector<Chain>::iterator it = chains.begin();it != chains.end();it++){
			const Edge& front = it->Edges.front();
			if(front.v1.OnTheEdge < 0){
				if(front.v2.OnTheEdge < 0){
					throw "front: All chain must starts from Edge";
				}
				else{
					std::swap(it->Edges.front().v1 , it->Edges.front().v2);
				}
			}

			const Edge& back = it->Edges.back();

			if(back.v2.OnTheEdge < 0){
				if(back.v1.OnTheEdge < 0){
					throw "back: All chain must starts from Edge";
				}
				else{
					std::swap(it->Edges.back().v1,it->Edges.back().v2);
				}
			}

			if(false && it->Edges.size() >=2 ){
				std::list<Edge>::iterator prev = it->Edges.begin();
				for(std::list<Edge>::iterator EdgeIt = it->Edges.begin()++; EdgeIt!=it->Edges.end() ;EdgeIt++,prev++){
					
					if(SameVertex(prev->v2.pointInTriangle , EdgeIt->v1.pointInTriangle)){
						continue;
					}else if(SameVertex(prev->v2.pointInTriangle , EdgeIt->v2.pointInTriangle)){
						std::swap(EdgeIt->v1,EdgeIt->v2);
					}
					else
						throw "Broken Chain";
				}
			}

		}
	}


	std::vector<Area> SortChains(CKContext* context){

		Chain lastChain[3];
		std::vector<Chain> chainsOnSameEdge[3];
		std::vector<Chain> EnclosedChains;

		std::vector<Area> ResultAreas;

		for(int i = 0;i < 3;i++){
			bool first = true;
			bool visible = triangle.visible[i];

			std::sort(chains.begin(),chains.end(),SortByGreater(i));

			//front(v1,v2) >>>>>> back(v1,v2)

			lastChain[i].valid = false;

			for(std::vector<Chain>::iterator it = chains.begin();
				it != chains.end();it++){
				if(it->type == ChainsOppositeToVertex && additionMsg2Vertex[it->additionalMsg] == i){
					//context->OutputToConsoleEx("ChainsOppositeToVertex");
					//context->OutputToConsoleEx("max[%d] = %f",i,it->max[i]);
					if(true == visible && true == first){
						//context->OutputToConsoleEx("true == visible && true == first");
						Area area;
						area.ParentNorm = triangle.norm;
						Edge E;
						

						VxVector& v1 = it->Edges.front().v1.point;
						VxVector& v2 = it->Edges.back().v2.point;

						VxVector norm = VxVectorCrossProduct(v2 - triangle.v[i], v1 - triangle.v[i]);
						norm.Normalize();

						
						//check direction

						if(!SameVertex(norm,triangle.norm)){
							it->reverse();
							
							std::list<Edge>::iterator end = it->Edges.end();
							for(std::list<Edge>::iterator Edgeit = it->Edges.begin(); Edgeit != it->Edges.end() ;Edgeit++){
								std::swap(Edgeit->v1,Edgeit->v2);
							}
							

							//context->OutputToConsoleEx("reject by norm ");
						    //context->OutputToConsoleEx("norm (%f,%f,%f) ParentNorm (%f,%f,%f)",norm.x,norm.y,norm.z,triangle.norm.x,triangle.norm.y,triangle.norm.z);

						}


						area.Boundary = *it;//just copy

						int startEdge = it->Edges.front().v1.OnTheEdge;
						int endEdge = it->Edges.back().v2.OnTheEdge;

						if(endEdge < 0 || startEdge < 0){
							int i = -1;
						}

						PointInTriangle _lastPoint = it->Edges.back().v2;

						for(std::vector<Chain>::reverse_iterator rit = chainsOnSameEdge[endEdge].rbegin(); rit != chainsOnSameEdge[endEdge].rend(); rit++ ){
							bool shouldReverse = false;
							if(rit->Edges.front().v1.pointInTriangle[i] > rit->Edges.back().v2.pointInTriangle[i]){
								rit->reverse();
								shouldReverse = true;
							}
							
							E.v1 = _lastPoint;
							if(shouldReverse)
								E.v2 = rit->Edges.front().v2;
							else
							    E.v2 = rit->Edges.front().v1;

							area.Boundary.Edges.push_back(E);

							for(std::list<Edge>::iterator _e = rit->Edges.begin(); _e != rit->Edges.end(); _e++ ){
								if(shouldReverse)
						            std::swap(_e->v1,_e->v2);
								area.Boundary.Edges.push_back(*_e);
							}

							_lastPoint = rit->Edges.back().v2;
						}

						chainsOnSameEdge[endEdge].clear();

						E.v1 = _lastPoint;

						E.v2 = PointInTriangle(triangle,triangle.v[i]);

						area.Boundary.Edges.push_back(E);

						_lastPoint = E.v2;


						for(std::vector<Chain>::iterator _it = chainsOnSameEdge[startEdge].begin(); _it != chainsOnSameEdge[startEdge].end();_it++ ){
							bool shouldReverse = false;
							if(_it->Edges.back().v2.pointInTriangle[i] > _it->Edges.front().v1.pointInTriangle[i]){
								_it->reverse();
								shouldReverse = true;
							}
							
							E.v1 = _lastPoint;
							if(shouldReverse)
								E.v2 = _it->Edges.front().v2;
							else
							    E.v2 = _it->Edges.front().v1;

							area.Boundary.Edges.push_back(E);

							for(std::list<Edge>::iterator _e = _it->Edges.begin(); _e != _it->Edges.end();_e++){
								if(shouldReverse)
						            std::swap(_e->v1,_e->v2);
								area.Boundary.Edges.push_back(*_e);
							}

							_lastPoint = _it->Edges.back().v2;
						}

						chainsOnSameEdge[startEdge].clear();

						E.v1 = _lastPoint;
						E.v2 = it->Edges.front().v1;

						area.Boundary.Edges.push_back(E);

						//context->OutputToConsoleEx("i %d",i);

						for(std::list<Edge>::iterator EdgeIt = area.Boundary.Edges.begin(); EdgeIt != area.Boundary.Edges.end();EdgeIt++){
							ShowVxVector(context,EdgeIt->v1.pointInTriangle);
					        ShowVxVector(context,EdgeIt->v2.pointInTriangle);
						}


						ResultAreas.push_back(area);
					}else if(false == visible && true == first){
						;
					}else if(true == visible && false == first){
						//DEBUGBREAK
						Area area;
						area.ParentNorm = triangle.norm;
						Edge E;

						VxVector& v1 = it->Edges.front().v1.point;
						VxVector& v2 = it->Edges.back().v2.point;

						VxVector norm = VxVectorCrossProduct(v2 - triangle.v[i], v1 - triangle.v[i]);
						norm.Normalize();

						
						//check direction

						if(!SameVertex(norm,triangle.norm)){
							it->reverse();
							
							std::list<Edge>::iterator end = it->Edges.end();
							for(std::list<Edge>::iterator Edgeit = it->Edges.begin(); Edgeit != it->Edges.end() ;Edgeit++){
								std::swap(Edgeit->v1,Edgeit->v2);
							}
							

							//context->OutputToConsoleEx("reject by norm ");
						    //context->OutputToConsoleEx("norm (%f,%f,%f) ParentNorm (%f,%f,%f)",norm.x,norm.y,norm.z,triangle.norm.x,triangle.norm.y,triangle.norm.z);

						}


						VxVector& lastv1 = lastChain[i].Edges.front().v1.point;
						VxVector& lastv2 = lastChain[i].Edges.back().v2.point;

						VxVector lastnorm = VxVectorCrossProduct(lastv2 - triangle.v[i], lastv1 - triangle.v[i]);
						lastnorm.Normalize();

						
						//check direction

						if(SameVertex(lastnorm,triangle.norm)){
							lastChain[i].Edges.reverse();
							
							std::list<Edge>::iterator end = lastChain[i].Edges.end();
							for(std::list<Edge>::iterator Edgeit = lastChain[i].Edges.begin(); Edgeit != lastChain[i].Edges.end() ;Edgeit++){
								std::swap(Edgeit->v1,Edgeit->v2);
							}
							

							//context->OutputToConsoleEx("reject by norm, direction should be inverse");
						    //context->OutputToConsoleEx("norm (%f,%f,%f) ParentNorm (%f,%f,%f)",lastnorm.x,lastnorm.y,lastnorm.z,triangle.norm.x,triangle.norm.y,triangle.norm.z);

						}


						
						Chain& nearChain = lastChain[i];
						Chain& farChain = *it;

						int startEdge = farChain.Edges.front().v1.OnTheEdge;
						int endEdge = farChain.Edges.back().v2.OnTheEdge;

						if(endEdge < 0 || startEdge < 0){
							int i = -1;
						}

						area.Boundary = farChain;//just copy

						PointInTriangle _lastPoint = farChain.Edges.back().v2;

						for(std::vector<Chain>::reverse_iterator rit = chainsOnSameEdge[endEdge].rbegin(); rit != chainsOnSameEdge[endEdge].rend(); rit++ ){
							bool shouldReverse = false;
							if(rit->Edges.front().v1.pointInTriangle[i] > rit->Edges.back().v2.pointInTriangle[i]){
								rit->reverse();
								shouldReverse = true;
							}
							
							E.v1 = _lastPoint;
							if(shouldReverse)
								E.v2 = rit->Edges.front().v2;
							else
							    E.v2 = rit->Edges.front().v1;

							area.Boundary.Edges.push_back(E);

							for(std::list<Edge>::iterator _e = rit->Edges.begin(); _e != rit->Edges.end(); _e++ ){
								if(shouldReverse)
									std::swap(_e->v1,_e->v2);
								area.Boundary.Edges.push_back(*_e);
							}

							_lastPoint = rit->Edges.back().v2;
						}

						chainsOnSameEdge[endEdge].clear();

						
						bool shouldReverse = false;
						if(nearChain.Edges.front().v1.OnTheEdge != endEdge){
							nearChain.reverse();
							shouldReverse = true;
						}

						E.v1 = _lastPoint;

						if(shouldReverse)
							E.v2 = nearChain.Edges.front().v2;
						else
						    E.v2 = nearChain.Edges.front().v1;

						area.Boundary.Edges.push_back(E);

						for(std::list<Edge>::iterator nit = nearChain.Edges.begin();nit != nearChain.Edges.end();nit++){
							if(shouldReverse)
						        std::swap(nit->v1,nit->v2);
							area.Boundary.Edges.push_back(*nit);
						}

						_lastPoint = nearChain.Edges.back().v2;



						for(std::vector<Chain>::iterator _it = chainsOnSameEdge[startEdge].begin(); _it != chainsOnSameEdge[startEdge].end();_it++ ){
							bool shouldReverse = false;
							if(_it->Edges.back().v2.pointInTriangle[i] > _it->Edges.front().v1.pointInTriangle[i]){
								_it->reverse();
								shouldReverse = true;
							}
							
							E.v1 = _lastPoint;
							if(shouldReverse)
								E.v2 = _it->Edges.front().v2;
							else
							    E.v2 = _it->Edges.front().v1;

							area.Boundary.Edges.push_back(E);

							for(std::list<Edge>::iterator _e = _it->Edges.begin(); _e != _it->Edges.end();_e++){
								if(shouldReverse)
						            std::swap(_e->v1,_e->v2);
								area.Boundary.Edges.push_back(*_e);
							}

							_lastPoint = _it->Edges.back().v2;
						}

						chainsOnSameEdge[startEdge].clear();

						E.v1 = _lastPoint;
						E.v2 = it->Edges.front().v1;

						area.Boundary.Edges.push_back(E);

						ResultAreas.push_back(area);
					}else if(false == visible && false == first){
						;
					}
					
					lastChain[i] = *it;
					lastChain[i].visible = visible;
					visible = !visible;

					//context->OutputToConsoleEx("lastChain[%d]->max[%d] = %f",i,i,lastChain[i].max[i]);

					first = false;
				}else if(it->type == ChainsOnSameEdge && it->additionalMsg == (i + 1)%3){
					chainsOnSameEdge[it->additionalMsg].push_back(*it);
					//context->OutputToConsoleEx("OnSameEdge");
				}
			}//for it
		}//for i

		//if(lastChain[1].valid)
		//context->OutputToConsoleEx("1 lastChain[1]->max[1] = %f",lastChain[1].max[1]);


		bool middleVisible;
		bool hasLastChain = false;
		for(int i = 0;i < 3;i++){
			if(lastChain[i].valid){
				hasLastChain = true;
				middleVisible = !lastChain[i].visible;
				//context->OutputToConsoleEx(" ~ lastChain[%d]->max[%d] = %f",i,i,lastChain[i].max[i]);
			}
		}

		if(!hasLastChain){
			middleVisible = triangle.visible[0];
			//context->OutputToConsoleEx("no Last Chain");
		}

		if(middleVisible){
			Area area;
			area.ParentNorm = triangle.norm;
			PointInTriangle _lastPoint;
			Edge E;

			//context->OutputToConsoleEx("middle visible");

			if(lastChain[1].valid){
				//context->OutputToConsoleEx("1");

				area.Boundary = lastChain[1];
				if(area.Boundary.Edges.back().v2.OnTheEdge != 0)
					area.Boundary.reverse();

				for(std::list<Edge>::iterator it = area.Boundary.Edges.begin(); it != area.Boundary.Edges.end();it++){
					std::swap(it->v1,it->v2);
				}

				_lastPoint = area.Boundary.Edges.back().v2;
				
			}else{
				_lastPoint = PointInTriangle(triangle,triangle.v[1]);
			}


			for(std::vector<Chain>::reverse_iterator rit = chainsOnSameEdge[0].rbegin(); rit != chainsOnSameEdge[0].rend(); rit++ ){
				bool shouldReverse = false;
				if(rit->Edges.front().v1.pointInTriangle[2] > rit->Edges.back().v2.pointInTriangle[2]){
					rit->reverse();
					shouldReverse = true;
				}

				E.v1 = _lastPoint;
				if(shouldReverse)
					E.v2 = rit->Edges.front().v2;
				else
					E.v2 = rit->Edges.front().v1;


				area.Boundary.Edges.push_back(E);

				for(std::list<Edge>::iterator _e = rit->Edges.begin(); _e != rit->Edges.end(); _e++ ){
					if(shouldReverse)
						std::swap(_e->v1,_e->v2);
					area.Boundary.Edges.push_back(*_e);
				}

				_lastPoint = rit->Edges.back().v2;

			}

			if(lastChain[2].valid){
				bool shouldReverse = false;
				if(lastChain[2].Edges.back().v2.OnTheEdge != 1){
					lastChain[2].reverse();
					shouldReverse = true;
				}
					
				E.v1 = _lastPoint;
				if(shouldReverse)
					E.v2 = lastChain[2].Edges.front().v2;
				else
				    E.v2 = lastChain[2].Edges.front().v1;

				area.Boundary.Edges.push_back(E);

				//context->OutputToConsoleEx("pre2");

				ShowVxVector(context,E.v1.pointInTriangle);
			    ShowVxVector(context,E.v2.pointInTriangle);

				for(std::list<Edge>::iterator _e = lastChain[2].Edges.begin(); _e != lastChain[2].Edges.end(); _e++ ){
					if(shouldReverse)
						std::swap(_e->v1,_e->v2);
					area.Boundary.Edges.push_back(*_e);
					//context->OutputToConsoleEx("l2");
					ShowVxVector(context,_e->v1.pointInTriangle);
					ShowVxVector(context,_e->v2.pointInTriangle);
				}

				_lastPoint = lastChain[2].Edges.back().v2;
			}else{
				E.v1 = _lastPoint;
				E.v2 = PointInTriangle(triangle,triangle.v[2]);

				area.Boundary.Edges.push_back(E);

				//context->OutputToConsoleEx("e2");

				ShowVxVector(context,E.v1.pointInTriangle);
			    ShowVxVector(context,E.v2.pointInTriangle);

				_lastPoint = E.v2;
			}


			for(std::vector<Chain>::iterator _it = chainsOnSameEdge[1].begin(); _it != chainsOnSameEdge[1].end();_it++ ){
				bool shouldReverse = false;
			    if(_it->Edges.back().v2.pointInTriangle[2] > _it->Edges.front().v1.pointInTriangle[2]){
						_it->reverse();
						shouldReverse = true;
				}
				E.v1 = _lastPoint;

				if(shouldReverse)
					E.v2 = _it->Edges.front().v2;
				else
					E.v2 = _it->Edges.front().v1;
				area.Boundary.Edges.push_back(E);

				for(std::list<Edge>::iterator _e = _it->Edges.begin(); _e != _it->Edges.end();_e++){
					if(shouldReverse)
						std::swap(_e->v1,_e->v2);
					area.Boundary.Edges.push_back(*_e);
				}

				_lastPoint = _it->Edges.back().v2;
			}

			
			if(lastChain[0].valid){
				bool shouldReverse = false;
				if(lastChain[0].Edges.back().v2.OnTheEdge != 2){
					lastChain[0].reverse();
					shouldReverse = true;
				}

				E.v1 = _lastPoint;

				if(shouldReverse)
					E.v2 = lastChain[0].Edges.front().v2;
				else
				    E.v2 = lastChain[0].Edges.front().v1;

				area.Boundary.Edges.push_back(E);

				//context->OutputToConsoleEx("pre0");


				ShowVxVector(context,E.v1.pointInTriangle);
				ShowVxVector(context,E.v2.pointInTriangle);


				for(std::list<Edge>::iterator _e = lastChain[0].Edges.begin(); _e != lastChain[0].Edges.end(); _e++ ){
					if(shouldReverse)
						std::swap(_e->v1,_e->v2);
					area.Boundary.Edges.push_back(*_e);
					//context->OutputToConsoleEx("l0");
					ShowVxVector(context,_e->v1.pointInTriangle);
					ShowVxVector(context,_e->v2.pointInTriangle);
				}

				_lastPoint = lastChain[0].Edges.back().v2;
			}else{
				E.v1 = _lastPoint;
				E.v2 = PointInTriangle(triangle,triangle.v[0]);

				area.Boundary.Edges.push_back(E);

				//context->OutputToConsoleEx("e0");

				ShowVxVector(context,E.v1.pointInTriangle);
			    ShowVxVector(context,E.v2.pointInTriangle);

				_lastPoint = E.v2;
			}


			for(std::vector<Chain>::reverse_iterator rit = chainsOnSameEdge[2].rbegin(); rit != chainsOnSameEdge[2].rend(); rit++ ){
				bool shouldReverse = false;
				if(rit->Edges.front().v1.pointInTriangle[1] > rit->Edges.back().v2.pointInTriangle[1]){
					rit->reverse();
					shouldReverse = true;
				}

				//context->OutputToConsoleEx("on same 2");

				E.v1 = _lastPoint;
				E.v2 = rit->Edges.front().v1;
				area.Boundary.Edges.push_back(E);

				ShowVxVector(context,E.v1.pointInTriangle);
				ShowVxVector(context,E.v2.pointInTriangle);

				for(std::list<Edge>::iterator _e = rit->Edges.begin(); _e != rit->Edges.end(); _e++ ){
					if(shouldReverse)
						std::swap(_e->v1,_e->v2);
					area.Boundary.Edges.push_back(*_e);
					ShowVxVector(context,_e->v1.pointInTriangle);
				    ShowVxVector(context,_e->v2.pointInTriangle);
				}

				_lastPoint = rit->Edges.back().v2;

			}

			//context->OutputToConsoleEx("close chain");

			E.v1 = _lastPoint;
			E.v2 = area.Boundary.Edges.front().v1;
			area.Boundary.Edges.push_back(E);



			/*

			if(lastChain[1] != chains.end()){
				E.v1 = _lastPoint;
				E.v2 = lastChain[1]->Edges.front().v1;
				area.Boundary.Edges.push_back(E);

				context->OutputToConsoleEx("connect to begin 1");

				ShowVxVector(context,E.v1.pointInTriangle);
			    ShowVxVector(context,E.v2.pointInTriangle);

			}else{
				E.v1 = _lastPoint;
				E.v2 = PointInTriangle(triangle,triangle.v[1]);
				area.Boundary.Edges.push_back(E);

				context->OutputToConsoleEx("e1");

				ShowVxVector(context,E.v1.pointInTriangle);
			    ShowVxVector(context,E.v2.pointInTriangle);
			}
			*/
			ResultAreas.push_back(area);
		}
		else{
			//context->OutputToConsoleEx("middle invisible");
			for(int i = 0;i < 3;i++){
				for(std::vector<Chain>::iterator it = chainsOnSameEdge[i].begin(); it != chainsOnSameEdge[i].end(); it++){
					//context->OutputToConsoleEx("chainsOnSameEdge");
					/*
					if(2 == i){
						if(it->Edges.back().v2.pointInTriangle[1] > it->Edges.front().v1.pointInTriangle[1]){
							it->reverse();
							for(std::list<Edge>::iterator _e = it->Edges.begin();_e != it->Edges.end();_e++){
								std::swap(_e->v1,_e->v2);
							}
						}
					}

					if(1 == i){
						if(it->Edges.back().v2.pointInTriangle[0] > it->Edges.front().v1.pointInTriangle[0]){
							it->reverse();
							for(std::list<Edge>::iterator _e = it->Edges.begin();_e != it->Edges.end();_e++){
								std::swap(_e->v1,_e->v2);
							}
						}
					}

					if(0 == i){
						if(it->Edges.back().v2.pointInTriangle[2] > it->Edges.front().v1.pointInTriangle[2]){
							it->reverse();
							for(std::list<Edge>::iterator _e = it->Edges.begin();_e != it->Edges.end();_e++){
								std::swap(_e->v1,_e->v2);
							}
						}
					}
					*/

					if(it->Edges.back().v2.pointInTriangle[(i + 3 - 1) % 3] > it->Edges.front().v1.pointInTriangle[(i + 3 - 1) % 3]){
						it->reverse();
						for(std::list<Edge>::iterator _e = it->Edges.begin();_e != it->Edges.end();_e++){
							std::swap(_e->v1,_e->v2);
						}
					}
					
					Area area;
					area.ParentNorm = triangle.norm;
					area.Boundary = *it;
					Edge E;
					E.v1 = area.Boundary.Edges.back().v2;
					E.v2 = area.Boundary.Edges.front().v1;

					area.Boundary.Edges.push_back(E);

					ResultAreas.push_back(area);
				}
			}
		}


		return ResultAreas;
	}//SortChains()
};

#endif