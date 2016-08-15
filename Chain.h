#ifndef __CHAIN_H__
#define __CHAIN_H__

#include <vector>
#include <list>
#include "MyUtility.h"
#include "Edge.h"
#include "Triangle.h"

enum ChainTypes{
	ChainsOnSameEdge,
	ChainsOppositeToVertex,
	EnclosedChains
};

class Chain
{
public:
	Chain(){
		max[0] = 0.0f;
		max[1] = 0.0f;
		max[2] = 0.0f;
		min[0] = 1.0f;
		min[1] = 1.0f;
		min[2] = 1.0f;

		valid = true;
	}

	void reverse(){
		Edges.reverse();
		/*
		std::list<Edge>::iterator end = Edges.end();
		for(std::list<Edge>::iterator it = Edges.begin(); it != end ;it++){
			std::swap(it->v1,it->v2);
		}
		*/
	}

	std::list<Edge> Edges;
	bool valid;
	bool lastTimeUsedForVisible;
	bool visible;
	ChainTypes type;
	int additionalMsg;
	float max[3];
	float min[3];
};

class ChainNode{
public:
	ChainNode():removed(false){
	}

	ChainNode* next;
	ChainNode* prev;

	bool removed;
	PointInTriangle V;
};

#endif