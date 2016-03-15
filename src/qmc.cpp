
#include "InputNg.h"
#include "Engine/InputCheck.h"
#include "Engine/Geometry.h"
#include "Models/HubbardOneOrbital/HubbardOneOrbital.h"
#include "Engine/Engine.h"
#include "Provenance.h"

const PsimagLite::String license=
"Copyright (c) 2015-2016, UT-Battelle, LLC\n"
"All rights reserved\n"
"\n"
"[Determinantal QMC++, Version 0.]\n"
"\n"
"---------------------------------------------------------\n"
"THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND\n"
"CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED\n"
"WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n"
"WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A\n"
"PARTICULAR PURPOSE ARE DISCLAIMED.\n"
"\n"
"Please see full open source license included in file LICENSE.\n"
"---------------------------------------------------------\n"
"\n";

typedef PsimagLite::InputNg<Qmc::InputCheck> InputNgType;
typedef Qmc::Geometry<InputNgType::Readable,double> GeometryType;
typedef Qmc::Model<GeometryType> ModelType;
typedef Qmc::Engine<ModelType> EngineType;

int main(int argc, char* argv[])
{
	Qmc::InputCheck inputCheck;
	PsimagLite::String filename="";
	int opt = 0;
	PsimagLite::String strUsage(argv[0]);
	strUsage += " -f filename";
	bool versionOnly = false;

	while ((opt = getopt(argc, argv,"f:V")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'V':
			versionOnly = true;
			break;
		default:
			inputCheck.usageMain(strUsage);
			return 1;
		}
	}

	bool b = (filename == "");
	if (b && !versionOnly) {
		std::cerr<<"USAGE is "<<argv[0]<<" -f filename | -V\n";
		return 1;
	}

	size_t npthreads = 1;
	PsimagLite::Concurrency concurrency(&argc,&argv,npthreads);

	// print license
	if (PsimagLite::Concurrency::root()) {
		std::cout<<license;
		Provenance provenance;
		std::cout<<provenance;
	}

	if (versionOnly) return 0;

	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	GeometryType geometry(io);
	ModelType model(io,geometry);

	EngineType engine(model,io);

	engine.main();
}
