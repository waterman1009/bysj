#ifndef _MLS_H
#define _MLS_H

#include "../Base/Mesh.h"
#include "../Object/Shape.h"
#include "../Base/Polygon.h"


class CMLS
{
	void setImage(CLattice* srcGrid,varray<CPloygon *> srcFeatureLine);

	void deform(varray<CPloygon *> targetFeatureLine);

	varray<Vec2> getTargetGrid();

	CLattice * _srcGrid;
	varray<CPloygon *> _srcFeatureLine;

	varray<Vec2> _targetGrid;

	varray<double> _preCalculatedArray;

	varray< varray<double> > _sumDeltaArray;

};

#endif