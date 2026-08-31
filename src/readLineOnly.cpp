#include "../include/header.h"
void readLineOnly(void) {
	//initialize variables
	int i;
	std::vector<double> x(0);
	std::vector<double> y(0);
	std::vector<double> z(0);
	std::vector<string> atomType;
	int atomNum = 0;
	string inputLineAtomInfo;
	std::vector<double> coreSite(4);

	//User input section
	cout << "Enter number of atoms: ";
	cin >> atomNum;
	cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	// Resize vectors to hold atom information
	x.resize(atomNum + 1);
	y.resize(atomNum + 1);
	z.resize(atomNum + 1);
	atomType.resize(atomNum + 1);

	// Loop to get atom information
	for (i = 1; i <= atomNum; i++) {
		cout << "Please Enter Atom " << i << "Info;Just like C " << endl;
		getline(cin, inputLineAtomInfo);
		AtomInfo atomInfos;
		if (extractSingleLine(inputLineAtomInfo, atomInfos)) {
			cout << "Atom Type: " << atomInfos.type << ", Coordinates: (" << atomInfos.x << ", " << atomInfos.y << ", " << atomInfos.z << ")" << endl;
			x[i] = atomInfos.x;
			y[i] = atomInfos.y;
			z[i] = atomInfos.z;
			atomType[i] = atomInfos.type;

		}
		else {
			cerr << "Error: Invalid atom information format." << endl;
			break;
		}

	}

	// Calculate core site
	coreSite = rWeightedSiteCalc(atomNum, atomType, x, y, z);
	cout << "Core Site Coordinates: (" << coreSite[0] << ", " << coreSite[1] << ", " << coreSite[2] << "," << coreSite[3] << ")" << endl;
	cout << "Welcome to the AOTM Core Site Next!" << endl;



}