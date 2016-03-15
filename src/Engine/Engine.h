#ifndef ENGINE_H
#define ENGINE_H
#include "EngineParams.h"
#include "Matrix.h"
#include "RandomForTests.h"

namespace Qmc {

template<typename ModelType>
class Engine {

	typedef typename ModelType::GeometryType::RealType RealType;
	typedef typename ModelType::IoInputType IoInputType;
	typedef PsimagLite::Matrix<RealType> MatrixType;
	typedef EngineParams<RealType,IoInputType> EngineParamsType;
	typedef typename PsimagLite::RandomForTests<RealType> RngType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorType;
	typedef std::pair<RealType,RealType> PairRealType;

	enum SpinEnum {SPIN_UP, SPIN_DOWN};

public:

	Engine(const ModelType& model,IoInputType& io)
	    : params_(io),
	      rng_(0),
	      model_(model),
	      X_(5),
	      eigs_(2),
	      spins_(params_.ntimes,model.geometry().numberOfSites())
	{
		spins_.setTo(1);

		SizeType nsites= spins_.n_col();
		for (SizeType i = 0; i < X_.size(); ++i)
			X_[i].resize(nsites,nsites);

		for (SizeType i = 0; i < eigs_.size(); ++i)
			eigs_[i].resize(nsites);

		expH0_ = model.geometry().matrix();
		std::cout<<expH0_;
		expH0_ *= (params_.beta/params_.ntimes);
		exp(expH0_);
		RealType alpha = params_.beta * model_.params().U / params_.ntimes;
		RealType arctan = atan(sqrt(abs(alpha)/4.0));

		if (alpha >= 0) {
			aUp_ = -alpha/2 - 2 * arctan;
			aDown_ = -alpha/2 + 2 * arctan;
		} else {
			aUp_ = -alpha/2 - 2 * arctan;
			aDown_ = -alpha/2 - 2 * arctan;
		}
	}

	void main()
	{
		for (SizeType i = 0; i < params_.thermalizations; ++i) {
			evolve(i);
			printSpins();
		}
	}

private:

	void evolve(SizeType iter)
	{
		PairRealType oldWeight = weightOf();
		SizeType nsites = model_.geometry().numberOfSites();
		RealType acceptance = 0;
		RealType density = 0;
		RealType sign = 0;

		for (SizeType time = 0; time < params_.ntimes; ++time) {
			for (SizeType site = 0; site < nsites; ++site) {
				flipSpin(time,site);
				calcX();
				PairRealType newWeight = weightOf();
				RealType ratio1 = newWeight.first/oldWeight.first;
				RealType ratio2 = newWeight.second/oldWeight.second;
//				std::cout<<"Old Weight= "<<oldWeight;
//				std::cout<<" New Weight= "<<newWeight<<" ratio= "<<ratio<<"\n";
				bool accept = acceptOrReject(fabs(ratio1*ratio2));
				int thisSign = (ratio1*ratio2 > 0) ? 1 : -1;
				if (accept) {
					oldWeight = newWeight;
					acceptance++;
					X_[2] = X_[3];
					eigs_[0] = eigs_[1];
				} else {
					flipSpin(time,site);
				}

				density += densityFunction()*thisSign;
				sign += thisSign;
			}
		}

		RealType factor = 1.0/(params_.ntimes * nsites);
		acceptance *= factor;
		std::cout<<"Acceptance ratio "<<acceptance<<"\n";

		sign *= factor;
		std::cout<<"Sign= "<<sign<<"\n";

		density *= factor;
		std::cout<<"Density= "<<(density/sign)<<"\n";
	}

	void flipSpin(SizeType time,SizeType site)
	{
		spins_(time,site) *= (-1);
	}

	PairRealType weightOf()
	{
		RealType ret1 = weightOf(X_[0]);
		RealType ret2 = weightOf(X_[1]);
//		int b1 = (ret1 > 0) ? 1 : -1;
//		int b2 = (ret2 > 0) ? 1 : -1;
//		if (b1 * b2 < 0) {
//			std::cout<<"Sign Problem!\n";
//		}

		return PairRealType(ret1,ret2);
	}

	RealType weightOf(const MatrixType& X)
	{
		X_[3] = X;
		diag(X_[3],eigs_[1],'V');

		RealType prod = 1.0;
		for (SizeType l = 0; l < eigs_[1].size(); ++l)
			prod *= (1.0 + eigs_[1][l]);

		return prod;
	}

	void calcX()
	{
		calcX(X_[0],SPIN_UP);
		calcX(X_[1],SPIN_DOWN);
	}

	void calcX(MatrixType& m,SpinEnum spin)
	{
		setToIdentity(m);
		MatrixType Xprev;
		SizeType nsites = model_.geometry().numberOfSites();
		VectorType expv(nsites);

		for (SizeType time = 0; time < params_.ntimes; ++time) {
			Xprev = m * expH0_;
			calcExpV(expv,time,spin);
			for (SizeType i = 0; i < nsites; ++i)
				for (SizeType j = 0; j < nsites; ++j)
					m(i,j) = Xprev(i,j) * expv[j];
		}
	}

	void calcExpV(VectorType& v,SizeType time,SpinEnum spin) const
	{
		RealType a = (spin == SPIN_UP) ? aUp_ : aDown_;
		for (SizeType i = 0; i < v.size(); ++i)
			v[i] = exp(a*spins_(time,i)-params_.mu*params_.beta/params_.ntimes);
	}

	bool acceptOrReject(RealType ratio) const
	{
		RealType prob = ratio/(1.0 + ratio);

		RealType r = rng_();
//		std::cout<<"prob = "<<prob<<" r= "<<r<<"\n";
		return (r<prob);
	}

	void setToIdentity(MatrixType& m) const
	{
		SizeType n = m.n_row();
		for (SizeType i = 0; i < n; ++i)
			for (SizeType j = 0; j < n; ++j)
				m(i,j) = (i == j) ? 1 : 0;
	}

	void printSpins() const
	{
		for (SizeType i = 0; i < spins_.n_row(); ++i) {
			for (SizeType j = 0; j < spins_.n_col(); ++j) {
				char c = (spins_(i,j) > 0) ? '+' : '-';
				std::cout<<c;
			}
			std::cout<<" ";
		}

		std::cout<<"\n";
	}

	RealType densityFunction()
	{
		computeGreenFunction();
		return trace(X_[4])/model_.geometry().numberOfSites();
	}

	RealType trace(const MatrixType& m) const
	{
		RealType sum = 0;
		SizeType n = m.n_row();
		for (SizeType i = 0; i < n; ++i)
			sum += m(i,i);

		return sum;
	}

	void computeGreenFunction()
	{
		SizeType n = X_[4].n_row();
		for (SizeType i = 0; i < n; ++i) {
			for (SizeType j = 0; j < n; ++j) {
				RealType sum = 0;
				for (SizeType k = 0; k < n; ++k) {
					sum += std::conj(X_[2](i,k)) * X_[2](j,k) / (1.0 + eigs_[0][k]);
				}

				X_[4](i,j) = -sum;
				if (i == j) X_[4](i,j) += 1.0;
			}
		}
	}

	EngineParamsType params_;
	mutable RngType rng_;
	const ModelType& model_;
	typename PsimagLite::Vector<MatrixType>::Type X_;
	typename PsimagLite::Vector<typename PsimagLite::Vector<RealType>::Type>::Type eigs_;
	PsimagLite::Matrix<int> spins_;
	RealType aUp_,aDown_;
	MatrixType expH0_;
}; // class Engine

} // namespace Qmc

#endif // ENGINE_H

