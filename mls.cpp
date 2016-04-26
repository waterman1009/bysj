#include "mls.h"
#include <cmath>



CMLS::CMLS() {}

CMLS::~CMLS() {}

void CMLS::setImage(CLattice* srcGrid,varray<CPloygon *> srcFeatureLine) {
	//	循环原图所有网格点
	_srcGrid = srcGrid;

	varray<double> arrayVertexPos = srcGrid->GetVertexPos();

	for (int i=0; i<arrayVertexPos.size()/2; k++) {
		Vec2 v(arrayVertexPos[ i+i ],arrayVertexPos[ i+i+1 ]);

		bool isBreak = false;

		Vec2 weightP;
		double sumP = 0.0;

		varray<double> deltaArray;
		//	对于每个点v，循环每条特征线，预计算参数
		for (int k=0; k<srcFeatureLine.size(); k++) {
			CPolygon* pPolygon = srcFeatureLine[k];
			int iVNumber = pPolygon->GetVNumber();
			bool isClosed = pPolygon->GetClosed();

			varray<Vec2>& arrayVertices = pPolygon->GetPolygon();

			int nEdge = isClosed ? iVNumber : iVNumber-1;
			int iEdge = 0;
			for(; iEdge < nEdge; iEdge ++) {
				int iV1 = iEdge, iV2 = (iV1+1)%iVNumber;

				//	得到line segments的两个端点Ai Bi
				Vec2 ai = arrayVertices[iV1],
					bi = arrayVertices[iV2];
				double deltaI = Cross( ai-v , ai-bi );
				double thltaI = atan(  ( Dot( bi-v , bi-ai ) )/( Cross( bi-v , bi-ai ) ) ) 
								- atan( ( Dot( ai-v , ai-bi ) )/( Cross( ai-v , ai-bi )) );
				double belta00 = Dot( ai-v , ai-v);
				double belta01 = Dot( ai-v , v-bi);
				double belta11 = Dot( v-bi , v-bi);

				double delta00,delta01,delta11;

				if ( fabs(deltaI)<1e-8 ) {
					// when deltaI ==0 
					if ( (ai._x<=v._x && v._x<=bi._x) || (bi._x<=v._x && v._x<=ai._x) ) {
						// when v is on the line segment defined by ai and bi
						isBreak = true;
						_precalculatedArray.push_back(ai._x);
						_precalculatedArray.push_back(ai._y);
						_precalculatedArray.push_back(bi._x);
						_precalculatedArray.push_back(bi._y);
					}
					else {
						// v is on the extension of one of these line segments
						delta00 = pow( (ai-bi).Magnitude(), 5.0 ) / 
							( 3.0 * pow( Dot( v-bi , bi-ai ), 1.0 ) * pow( Dot( ai-v , bi-ai ), 3.0 ) );
						delta01 = -pow( (ai-bi).Magnitude(), 5.0 ) /
							( 6.0 * pow( Dot( v-bi , bi-ai ), 2.0 ) * pow( Dot( ai-v , bi-ai ), 2.0 ) );
						delta11 = pow( (ai-bi).Magnitude(), 5.0 ) /
							( 3.0 * pow( Dot( v-bi , bi-ai ), 3.0 ) * pow( Dot( ai-v , bi-ai ), 1.0 ) );
					}
				}
				else {
					delta00 = (ai-bi).Magnitude() / (2.0 * deltaI * deltaI) * 
						(belta01/belta00 - belta11*thltaI/deltaI );
					delta01 = (ai-bi).Magnitude() / (2.0 * deltaI * deltaI) *
						(1.0 - belta01*thltaI/deltaI);
					delta11 = (ai-bi).Magnitude() / (2.0 * deltaI * deltaI) *
						(belta01/belta11 - belta00*thltaI/deltaI );
				}

				weightP += ai*( delta00 + delta01 ) + bi*( delta01 + delta11 );
				sumP += delta00 + delta01 + delta01 + delta11;
				deltaArray.push_back(isBreak==true?1.0:0.0);
				deltaArray.push_back(delta00);
				deltaArray.push_back(delta01);
				deltaArray.push_back(delta11);
			}

			if (isBreak) {
				break;
			}
		}

		_sumDeltaArray.push_back(deltaArray);

		if (isBreak) {
			continue;
		}

		int deltaPos = 0;

		weightP /= sumP;

		_precalculatedArray.push_back(weightP._x);
		_precalculatedArray.push_back(weightP._y);

		// isBreak is false

		for (int k=0; k<srcFeatureLine.size(); k++) {
			CPolygon* pPolygon = srcFeatureLine[k];
			int iVNumber = pPolygon->GetVNumber();
			bool isClosed = pPolygon->GetClosed();

			varray<Vec2>& arrayVertices = pPolygon->GetPolygon();

			int nEdge = isClosed ? iVNumber : iVNumber-1;
			int iEdge = 0;
			for(; iEdge < nEdge; iEdge ++) {
				int iV1 = iEdge, iV2 = (iV1+1)%iVNumber;
				//	得到line segments的两个端点Ai Bi
				Vec2 ai = arrayVertices[iV1] - weightP,
					bi = arrayVertices[iV2] - weightP;
				deltaPos++;
				double delta00 = deltaArray[deltaPos++];
				double delta01 = deltaArray[deltaPos++];
				double delta11 = deltaArray[deltaPos++];

				double a00 = delta00 * ai._x + delta01 * bi._x;
				double a01 = delta00 * ai._y + delta01 * bi._y;
				double a10 = delta00 * ai._y + delta01 * bi._y;
				double a11 = delta00 * -ai._x + delta01 * -bi._x;
				double a20 = delta01 * ai._x + delta11 * bi._x;
				double a21 = delta01 * ai._y + delta11 * bi._y;
				double a30 = delta01 * ai._y + delta11 * bi._y;
				double a31 = delta01 * -ai._x + delta11 * -bi._x;

				Vec2 vp = v - weightP;
				
				_precalculatedArray.push_back( a00*vp._x + a01*vp._y );
				_precalculatedArray.push_back( a00*vp._y + a01*-vp._x );
				_precalculatedArray.push_back( a10*vp._x + a11*vp._y );
				_precalculatedArray.push_back( a10*vp._y + a11*-vp._x );
				_precalculatedArray.push_back( a20*vp._x + a21*vp._y );
				_precalculatedArray.push_back( a20*vp._y + a21*-vp._x);
				_precalculatedArray.push_back( a30*vp._x + a31*vp._y);
				_precalculatedArray.push_back( a30*vp._y + a31*-vp._x);
			}
		}
	}

	return ;
}

void CMLS::deform(varray<CPloygon *> targetFeatureLine) {
	//	循环原图所有网格点
	varray<double> arrayVertexPos = _srcGrid->GetVertexPos();

	int precalucatedPos = 0;
	for (int i=0; i<arrayVertexPos.size()/2; i++) {
		Vec2 v(arrayVertexPos[ i+i ],arrayVertexPos[ i+i+1 ]);
		varray<double> deltaArray = _sumDeltaArray[i];

		int deltaPos = 0;
		bool isBreak = false;
		Vec2 weightQ;
		double sumQ = 0.0;

		for (int k=0; k<targetFeatureLine.size(); k++) {
			CPolygon* pPolygon = targetFeatureLine[k];
			int iVNumber = pPolygon->GetVNumber();
			bool isClosed = pPolygon->GetClosed();

			varray<Vec2>& arrayVertices = pPolygon->GetPolygon();

			int nEdge = isClosed ? iVNumber : iVNumber-1;
			int iEdge = 0;
			for(; iEdge < nEdge; iEdge ++) {
				int iV1 = iEdge, iV2 = (iV1+1)%iVNumber;

				//	得到line segments的两个端点Ci Di
				Vec2 ci = arrayVertices[iV1],
					di = arrayVertices[iV2];

				if ( fabs(deltaArray[deltaPos++])>1e-8 ) {
					// v is on the line segments defined by ai and bi
					double ax = _precalculatedArray[precalucatedPos++];
					double ay = _precalculatedArray[precalucatedPos++];
					double bx = _precalculatedArray[precalucatedPos++];
					double by = _precalculatedArray[precalucatedPos++];
					Vec2 ai(ax,ay);
					Vec2 bi(bx,by);
					double fvx = (di._x-ci._x)*(v._x-ai._x)/(bi._x-ai._x)+ci._x;
					double fvy = (di._y-ci._y)*(v._y-ai._y)/(bi._y-ai._y)+ci._y;
					_targetGrid.push_back(new Vec2(fvx,fvy));
					break;
				}

				double delta00 = deltaArray[deltaPos++];
				double delta01 = deltaArray[deltaPos++];
				double delta11 = deltaArray[deltaPos++];

				weightQ += ci*(delta00+delta01)+di*(delta01+delta11);
				sumQ += delta00+delta01+delta01+delta11;

			}
			if (isBreak) {
				break;
			}
		}

		if (isBreak) {
			continue;
		}

		weightQ /= sumQ;

		double px = _precalculatedArray[precalucatedPos++];
		double py = _precalculatedArray[precalucatedPos++];
		Vec2 weightP(px,py);

		Vec2 frv;
		
		for (int k=0; k<targetFeatureLine.size(); k++) {
			CPolygon* pPolygon = targetFeatureLine[k];
			int iVNumber = pPolygon->GetVNumber();
			bool isClosed = pPolygon->GetClosed();

			varray<Vec2>& arrayVertices = pPolygon->GetPolygon();

			int nEdge = isClosed ? iVNumber : iVNumber-1;
			int iEdge = 0;
			for(; iEdge < nEdge; iEdge ++) {
				int iV1 = iEdge, iV2 = (iV1+1)%iVNumber;
				//	得到line segments的两个端点Ci Di
				Vec2 ci = arrayVertices[iV1]-weightQ,
					di = arrayVertices[iV2]-weightQ;
				double a00 = _precalculatedArray[precalucatedPos++];
				double a01 = _precalculatedArray[precalucatedPos++];
				double a10 = _precalculatedArray[precalucatedPos++];
				double a11 = _precalculatedArray[precalucatedPos++];
				double a20 = _precalculatedArray[precalucatedPos++];
				double a21 = _precalculatedArray[precalucatedPos++];
				double a30 = _precalculatedArray[precalucatedPos++];
				double a31 = _precalculatedArray[precalucatedPos++];

				frv += new Vec2(ci._x*a00+ci._y*a10+di._x*a20+di._y*a30 , ci._x*a01+ci._y*a11+di._x*a21+di._y*a31);

			}
		}
		Vec2 ans = frv/frv.Magnitude()*(v-weightP).Magnitude()+weightQ;

		_targetGrid.push_back(ans);
	}
	return ;
}

varray<Vec2> CMLS::getTargetGrid() {
	return _targetGrid;
}