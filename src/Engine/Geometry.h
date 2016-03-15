#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "Matrix.h"
#include "../../../FreeFermions/src/GeometryParameters.h"
#include "../../../FreeFermions/src/GeometryLibrary.h"

namespace Qmc {

template<typename IoInputType_,typename RealType_>
class Geometry {

public:

	typedef RealType_ RealType;
	typedef IoInputType_ IoInputType;
	typedef PsimagLite::Matrix<RealType> MatrixType;
	typedef FreeFermions::GeometryParameters<RealType,IoInputType> GeometryParamsType;

	Geometry(IoInputType& io)
	    : geometryParams_(io),data_(geometryParams_)
	{}

	SizeType numberOfSites() const
	{
		return data_.matrix().n_row();
	}

	const MatrixType& matrix() const
	{
		return data_.matrix();
	}

private:

	GeometryParamsType geometryParams_;
	FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> data_;
}; // class Geometry

}
#endif // GEOMETRY_H

