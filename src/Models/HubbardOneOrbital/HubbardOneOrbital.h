#ifndef HUBBARDONEORBITAL_H
#define HUBBARDONEORBITAL_H
#include "ParamsHubbardOneOrbital.h"

namespace Qmc {

template<typename GeometryType_>
class Model {

public:

	typedef GeometryType_ GeometryType;
	typedef typename GeometryType::IoInputType IoInputType;
	typedef typename GeometryType::RealType RealType;
	typedef ParamsHubbardOneOrbital<IoInputType,RealType>  ModelParamsType;

	Model(IoInputType& io,GeometryType& geometry)
	    : geometry_(geometry), params_(io)
	{}

	const GeometryType& geometry() const { return geometry_; }

	const ModelParamsType& params() const { return params_; }

private:

	const GeometryType& geometry_;
	ModelParamsType params_;
};

}

#endif // HUBBARDONEORBITAL_H
