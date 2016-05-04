#ifndef _MLS_H
#define _MLS_H

#include "../base/Mesh.h"
#include "../base/Vec2.h"
#include "../Base/Lattice.h"
#include "../Base/Polygon.h"
#include <math.h>


class CMLS
{
public:

	CMLS();
	~CMLS();

	void SetLatticeFeatureLine(CLattice* srcGrid,const varray<CPolygon *>& srcFeatureLine);

	void DeformMesh(const varray<CPolygon *>& targetFeatureLine);

	const varray<double>& GetDeformLatticeVPos() const;

private:

	CLattice * _srcGrid;
	varray<CPolygon *> _srcFeatureLine;

	varray<double> _targetGrid;

	varray<double> _precalculatedArray;

	varray< varray<double> > _sumDeltaArray;

};

#endif