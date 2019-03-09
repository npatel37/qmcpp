#ifndef ENGINE_H
#define ENGINE_H
#include "EngineParams.h"
#include "Matrix.h"
#include "RandomForTests.h"

/*
 * X_ = vector matrix of size 5 x nsites x nsites
 *      X_[0] - O spin up --> B_M B_M-1 .... B_1  with  B_l = e^{-dT K} e^{-dT V(l)}             det{O(snew)}/det{O(sold)}
 *      X_[1] - O spin down --> same as X_[0]
 *      X_[2] -
 *      X_[3] -
 *      X_[4] - Green's function
*/

/* need to find a way to store X_[0] and X_[1] for all time steps!
 * Needed for calculations of the green's function.
*/

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
          GreenUp_(params_.ntimes),
          GreenDn_(params_.ntimes),
          spins_(params_.ntimes,model.geometry().numberOfSites()) // -- HS field - nsites x #time-slice
	{
		spins_.setTo(1);

		SizeType nsites= spins_.n_col();
		for (SizeType i = 0; i < X_.size(); ++i)
			X_[i].resize(nsites,nsites);

		for (SizeType i = 0; i < eigs_.size(); ++i)
			eigs_[i].resize(nsites);

        for (SizeType time = 0; time < params_.ntimes; ++time) {
            GreenUp_[time].resize(nsites,nsites);
            GreenDn_[time].resize(nsites,nsites);
        }

        expH0_ = model.geometry().matrix();
		std::cout<<expH0_;
		expH0_ *= (params_.beta/params_.ntimes);
        exp(expH0_);                                    // --- non-interacting part of the propogator e^{-dT K}
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


        SizeType nsites = model_.geometry().numberOfSites();
        SizeType ntimes = params_.ntimes;
		for (SizeType i = 0; i < params_.thermalizations; ++i) {
			evolve(i);
            RealType outmu = chemicalpotential2(params_.mu,params_.filling);
            params_.mu = outmu;
			printSpins();
            //cout.flush();
		}



        std::cout << " === check: density as a function of time === " << std::endl;
        computeGreenFunctionAll();
        for (SizeType time=0; time<ntimes; time++) {
            RealType densityup = 0.0, densitydn = 0.0;
            for (SizeType i = 0; i < nsites; ++i) {
                densityup += GreenUp_[time](i,i);
                densitydn += GreenDn_[time](i,i);
            }
            std::cout << time << " " << densityup << " " << densitydn << std::endl;
        }

        //RealType outmu = chemicalpotential(0.0,1.0);
        //std::cout << outmu << std::endl;



        std::cout << " === check band-width: chemical potential sweep === " << std::endl;
        for(SizeType i=0; i<=100;i++) {
            params_.mu = -6.0 + i*0.12;
            computeGreenFunctionAll();
            for (SizeType time=ntimes-1; time<ntimes; time++) {
                RealType densityup = 0.0, densitydn = 0.0;
                for (SizeType i = 0; i < nsites; ++i) {
                    densityup += GreenUp_[time](i,i);
                    densitydn += GreenDn_[time](i,i);
                }
                std::cout << params_.mu << " " << densityup/nsites << " " << densitydn/nsites << std::endl;
            }

        }


	}

private:

    RealType chemicalpotential2(RealType muin,double filling){
        //double temp;
        SizeType nsites = model_.geometry().numberOfSites();
        //SizeType finaltime = params_.ntimes-1;
        RealType N=filling;
        RealType diff, dcen, mucen;
        RealType muTolerance = 1e-5;


        RealType mulow=-model_.params().U-4.0;
        params_.mu = mulow;
        computeGreenFunctionAll();
        RealType dlow = (densityup() + densitydown())/nsites;

        RealType muhigh=model_.params().U+4.0;
        params_.mu = muhigh;
        computeGreenFunctionAll();
        RealType dhigh = (densityup() + densitydown())/nsites;


        for(int i=0;i<200;i++){

            mucen = 0.5*(mulow+muhigh);
            params_.mu = mucen;
            computeGreenFunctionAll();
            RealType nup = densityup()/nsites;
            RealType ndown = densitydown()/nsites;
            dcen = nup + ndown;

            diff = fabs(N-dcen);

            if(diff<=muTolerance){
               break;
            }
            else if(dcen>N){
                muhigh=mucen;
                dhigh = dcen;
            }
            else if(dcen<N){
                mulow = mucen;
                dlow = dcen;
            }
        }

        muin = mucen;
        params_.mu = muin;

        std::cout << "mu adjusted = " << muin
                  << " Numb elec (up,down) = (" << densityup() << "," << densitydown() << ")"
                  << " total density = " << dcen/nsites
                  << std::endl;
        assert(diff<=muTolerance);
        return muin;
    } // ----------

    RealType chemicalpotential1(RealType muin,double filling){
        //double temp;
        SizeType nsites = model_.geometry().numberOfSites();
        RealType nup, ndown, N1, diff;
        RealType N=filling*nsites;
        muin = params_.mu;
        RealType muTolerance = 1e-5;

        for(int i=0;i<200;i++){             // do not iterate more than 100 times - too expensive!
            computeGreenFunctionAll();
            nup = densityup();
            ndown = densitydown();
            N1=nup+ndown;
            diff = fabs(N-N1);

            if(diff<=muTolerance){
                break;
            }
            else if(N1>N){
                //params_.mu = params_.mu*(1+diff*0.01);
                params_.mu = params_.mu*(1+diff*0.01);
                muin = params_.mu;
            }
            else if(N1<N){
                //params_.mu = params_.mu*(1-diff*0.01);
                params_.mu = params_.mu*(1-diff*0.01);
                muin = params_.mu;
            }
        }

        std::cout << "mu adjusted = " << muin
                  << " Numb elec (up,down) = (" << nup << "," << ndown << ")"
                  << " total density = " << N1/nsites
                  << std::endl;
        assert(diff<=muTolerance);
        return muin;
    } // ----------




    RealType chemicalpotentialOLD(RealType muin,double filling){
        //double temp;
        SizeType nsites = model_.geometry().numberOfSites();
        //SizeType finaltime = params_.ntimes-1;
        RealType N=filling;
        RealType nup, ndown, N1, diff;
        RealType mulow=-model_.params().U-4.0;
        RealType muhigh=model_.params().U+4.0;
        RealType muTolerance = 1e-5;

        for(int i=0;i<200;i++){

            computeGreenFunctionAll();
            nup = densityup()/nsites;
            ndown = densitydown()/nsites;
            N1=(nup+ndown);
            diff = fabs(N-N1);

            if(diff<=muTolerance){
               break;
            }
            else if(N1>N){
                muhigh=muin;
                muin=0.5*(mulow+muin);
                params_.mu = muin;
            }
            else if(N1<N){
                mulow=muin;
                muin=0.5*(muhigh+muin);
                params_.mu = muin;
            }
        }

        std::cout << "mu adjusted = " << muin
                  << " Numb elec (up,down) = (" << nup << "," << ndown << ")"
                  << " total density = " << N1/nsites
                  << std::endl;
        assert(diff<=muTolerance);
        return muin;
    } // ----------

    RealType densityup(){
        RealType out=0.0;
        SizeType nsites = model_.geometry().numberOfSites();
        SizeType finaltime = params_.ntimes-1;

        for (SizeType i = 0; i < nsites; ++i) {
            out += GreenUp_[finaltime](i,i);
        }

        return out;
    }

    RealType densitydown(){
        RealType out=0.0;
        SizeType nsites = model_.geometry().numberOfSites();
        SizeType finaltime = params_.ntimes-1;

        for (SizeType i = 0; i < nsites; ++i) {
            out += GreenDn_[finaltime](i,i);
        }

        return out;
    }

    void setGreenZero() {
        SizeType nsites = model_.geometry().numberOfSites();
        SizeType ntimes = params_.ntimes;
        for (SizeType time=0; time<ntimes; time++) {
            for (SizeType i = 0; i < nsites; ++i) {
                for (SizeType j = 0; j < nsites; ++j) {
                    GreenUp_[time](i,j) = 0.0;
                    GreenDn_[time](i,j) = 0.0;
                }
            }
        }
    }

    void computeGreenFunctionAll() {

        setGreenZero();
        for (SizeType time = 0; time < params_.ntimes; ++time) {

            SizeType tp1 = time+1;

            // -- spin up
            calcX(X_[0],SPIN_UP,tp1);
            diag(X_[0],eigs_[0],'V');

            // -- spin down
            calcX(X_[1],SPIN_DOWN,tp1);
            diag(X_[1],eigs_[1],'V');

            // -- calculate the Green's Function at time
            SizeType n = X_[0].n_row();
            for (SizeType i = 0; i < n; ++i) {
                for (SizeType j = 0; j < n; ++j) {
                    RealType sumup=0, sumdn=0;
                    for (SizeType k = 0; k < n; ++k) {
                        sumup += PsimagLite::conj(X_[0](i,k)) * X_[0](j,k) / (1.0 + eigs_[0][k]);
                        sumdn += PsimagLite::conj(X_[1](i,k)) * X_[1](j,k) / (1.0 + eigs_[1][k]);
                    }

                    GreenUp_[time](i,j) = -sumup;
                    if (i == j) GreenUp_[time](i,j) += 1.0;
                    GreenDn_[time](i,j) = -sumdn;
                    if (i == j) GreenDn_[time](i,j) += 1.0;
                }
            }

        }

    }

	void evolve(SizeType iter)
	{
        diag(X_[3],eigs_[1],'V');
		PairRealType oldWeight = weightOf();
		SizeType nsites = model_.geometry().numberOfSites();
		RealType acceptance = 0;
		RealType density = 0;
		RealType sign = 0;

		for (SizeType time = 0; time < params_.ntimes; ++time) {
			for (SizeType site = 0; site < nsites; ++site) {
				flipSpin(time,site);
                calcX(params_.ntimes);                                        // -- calculate X_[0] and X_[1]
                PairRealType newWeight = weightOf();            // -- determinant of spin up X[0] and spin dn X[1]
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
                    flipSpin(time,site);                        // -- if reject, flip the spin back to original conf.
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
        RealType ret1 = weightOf(X_[0]); // det[1+X0] = det[1 + Bm Bm-1 ... B1] spin up
        RealType ret2 = weightOf(X_[1]); // det[1+X1] = det[1 + Bm Bm-1 ... B1] spin dn
//		int b1 = (ret1 > 0) ? 1 : -1;
//		int b2 = (ret2 > 0) ? 1 : -1;
//		if (b1 * b2 < 0) {
//			std::cout<<"Sign Problem!\n";
//		}

		return PairRealType(ret1,ret2);
	}

    // -- determinant of matrix 1+X
    // -- der[diagonal] matrix is product of the diagonal
	RealType weightOf(const MatrixType& X)
	{
		X_[3] = X;
		diag(X_[3],eigs_[1],'V');

		RealType prod = 1.0;
		for (SizeType l = 0; l < eigs_[1].size(); ++l)
			prod *= (1.0 + eigs_[1][l]);

		return prod;
	}

    void calcX(SizeType& maxtime)
	{
        calcX(X_[0],SPIN_UP,maxtime);
        calcX(X_[1],SPIN_DOWN,maxtime);
	}

    // B_M B_M-1 .... B_1
    // B_l = e^{-dT K} e^{-dT V(l)}
    void calcX(MatrixType& m,SpinEnum spin,SizeType& maxtime)
	{
		setToIdentity(m);
		MatrixType Xprev;
		SizeType nsites = model_.geometry().numberOfSites();
		VectorType expv(nsites);

        for (SizeType time = 0; time < maxtime; ++time) {
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
            v[i] = exp(a*spins_(time,i)+params_.mu*params_.beta/params_.ntimes);
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
					sum += PsimagLite::conj(X_[2](i,k)) * X_[2](j,k) / (1.0 + eigs_[0][k]);
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

    typename PsimagLite::Vector<MatrixType>::Type GreenUp_;
    typename PsimagLite::Vector<MatrixType>::Type GreenDn_;

    PsimagLite::Matrix<int> spins_;
	RealType aUp_,aDown_;
	MatrixType expH0_;
}; // class Engine

} // namespace Qmc

#endif // ENGINE_H

