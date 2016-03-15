#ifndef PARAMSHUBBARDONEORBITAL_H
#define PARAMSHUBBARDONEORBITAL_H

namespace Qmc {

template<typename IoInputType,typename RealType>
struct ParamsHubbardOneOrbital {

	ParamsHubbardOneOrbital(IoInputType& io)
	{
		io.readline(U,"HubbardU=");
	}

	RealType U;
};

}
#endif // PARAMSHUBBARDONEORBITAL_H

