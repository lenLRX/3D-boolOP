#ifndef __RANGE_H__
#define __RANGE_H__

template<typename T>
class Range{
public:
	Range(T b1,T b2){
		if(b1 < b2){
			lb = b1;
			ub = b2;
		}else{
			ub = b1;
			lb = b2;
		}
	}

	T lb;//lower Bound
	T ub;//upper Bound
};

#endif