//////////////////////////////////////////////////////////////////////////////////////////////////////////
//		            MyVertexManager
//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "CKAll.h"
#include "MyVertexManager.h"

MyVertexManager::MyVertexManager(CKContext *Context) : 
CKBaseManager(Context, MyVertexManagerGUID, "MyVertexManager")
{
	Context->RegisterNewManager(this);
}

MyVertexManager::~MyVertexManager()
{
}

//---  Called at the beginning of each process loop.
CKERROR MyVertexManager::PreProcess()
{
	return CKERR_NOTIMPLEMENTED;
}

//---  Called at the end of each process loop.
CKERROR MyVertexManager::PostProcess()
{
	return CKERR_NOTIMPLEMENTED;
}




