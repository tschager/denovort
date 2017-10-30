# Improved De Novo Peptide Sequencing using LC Retention Time Information

## Installation

Download and extract the source code of OpenMS 2.0.1 (https://github.com/OpenMS/OpenMS/releases).
Copy the source code of this directory in the appropriate folders

	rsync -avhpI src OpenMS-Release2.0.1/src

Build OpenMS (see https://github.com/OpenMS/OpenMS/wiki/Building-OpenMS)

	mkdir OpenMS-build
	cd OpenMS-build
	cmake -DOPENMS_CONTRIB_LIBS="/PATH/TO/contrib-build"  -DBOOST_USE_STATIC=OFF ../OpenMS-Release2.0.1/

## Usage

Examples for the configuration files are in the folders parameters/ and coefs/.
The algorithm is used as follows

	./Openms-build/bin/DeNovoSymDiff -ini parameters/parameters_lin.ini -in spectra/ms2.dta2d -out test.idXML -parent_mass 1195.617064 -algorithm:annotation_sequence QAANQQTVEAK -algorithm:rho 0.9
