#pragma once
#include "shared.h"
#include "gridmath.h"

class MarchingCubes {
private:
	//(i,j,k) toward (i+1,j,k),(i,j+1,k),(i,j,k+1), vertex on edge
	aryi vertlisx, vertlisy, vertlisz;
public:
	vector<Float> verts;
	vector<int> faces;
	MarchingCubes(void);
	~MarchingCubes();
	void List_Verts(const aryf & f, Float isoval);
	void List_Faces(const aryf & f, Float isoval);
	void Reconstruct(const aryf &f, Float isoval);
	void Dump_Obj(const char *filename);
	void Dump_GOC(const char * gocname, const char *pngname, int outputW, int outputH);//GOC: scene file format used in wmd's renderer
};

