
#ifndef INPUT_CHECK_H
#define INPUT_CHECK_H
#include <vector>
#include <stdexcept>
#include "Options.h"

namespace Qmc {

	class InputCheck {

		typedef PsimagLite::Options::Readable OptionsReadableType;

	public:

		InputCheck() : optsReadable_(0) {}

		~InputCheck()
		{
			if (optsReadable_!=0) delete optsReadable_;
		}

		bool check(const PsimagLite::String& label,const PsimagLite::Vector<PsimagLite::String>::Type& vec,size_t line) const
		{
			if (label=="JMVALUES") {
				if (vec.size()!=2) return error1("JMVALUES",line);
				return true;
			} else if (label=="RAW_MATRIX") {
				size_t row = atoi(vec[0].c_str());
				size_t col = atoi(vec[1].c_str());
				size_t n = row*col;
				if (vec.size()!=n+2) return error1("RAW_MATRIX",line);
				return true;
			} else if (label=="Connectors") {
				return true;
			} else if (label=="MagneticField") {
				return true;
			} else if (label=="FiniteLoops") {
				size_t n = atoi(vec[0].c_str());
				if (vec.size()!=3*n+1)  return error1("FiniteLoops",line);
				return true;
			}
			return false;
		}

		void check(const PsimagLite::String& label,const PsimagLite::String& val,size_t line)
		{
			if (label!="SolverOptions") return;
			PsimagLite::Vector<PsimagLite::String>::Type registerOpts;

			registerOpts.push_back("none");

			PsimagLite::Options::Writeable optWriteable(registerOpts,PsimagLite::Options::Writeable::PERMISSIVE);
			optsReadable_ = new  OptionsReadableType(optWriteable,val);
		}

		bool isSet(const PsimagLite::String& thisOption) const
		{
			return optsReadable_->isSet(thisOption);
		}

		void checkForThreads(size_t nthreads) const
		{
			if (nthreads==1) return;

			PsimagLite::String message1(__FILE__);
			message1 += " FATAL: You are requesting nthreads>0 but you did not compile with USE_PTHREADS enabled\n";
			message1 += " Either set Threads=1 in the input file (you won't have threads though) or\n";
			message1 += " add -DUSE_PTHREADS to the CPP_FLAGS in your Makefile and recompile\n";
			throw PsimagLite::RuntimeError(message1.c_str());
		}

		void usageMain(const PsimagLite::String& name) const
		{
			std::cerr<<"USAGE is "<<name<<"\n";
		}

	private:

		bool error1(const PsimagLite::String& message,size_t line) const
		{
			PsimagLite::String s(__FILE__);
			s += " : Input error for label " + message + " near line " + ttos(line) + "\n";
			throw PsimagLite::RuntimeError(s.c_str());

		}

		OptionsReadableType* optsReadable_;

	}; // class InputCheck
} // namespace Qmc

/*@}*/
#endif

