#ifndef ENGINEPARAMS_H
#define ENGINEPARAMS_H

namespace Qmc {

template<typename RealType, typename IoInputType>
struct EngineParams {

	EngineParams(IoInputType& io)
	{
		io.readline(thermalizations,"Thermalizations=");
		io.readline(ntimes,"NumberOfTimes=");
		io.readline(beta,"Beta=");
		io.readline(mu,"ChemicalPotential=");
	}

	SizeType thermalizations;
	SizeType ntimes;
	RealType beta;
	RealType mu;
};

}
#endif // ENGINEPARAMS_H
