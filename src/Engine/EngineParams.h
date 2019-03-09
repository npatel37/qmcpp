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
        dtau = beta/ntimes;
		io.readline(mu,"ChemicalPotential=");
        io.readline(filling,"filling=");
	}

	SizeType thermalizations;
	SizeType ntimes;
	RealType beta;
    RealType dtau;
	RealType mu;
    RealType filling;
};

}
#endif // ENGINEPARAMS_H

