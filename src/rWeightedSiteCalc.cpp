#include "../include/header.h"
std::vector<double> rWeightedSiteCalc(int atomNum, std::vector<string> atomType, std::vector<double> x, std::vector<double> y, std::vector<double> z) {
	//Declare variables
	int i;
	double r = 0;
	double rCube = 0;
	double xTotal = 0;
	double yTotal = 0;
	double zTotal = 0;
	double rCubeTotal = 0;
	std::vector<double> coreSite(4);
	//Atom types could be "H", "C", "N", "O", "S" etc.

	if (atomNum > 0 && initArrayDect(atomNum, atomType, x, y, z)) {
		for (i = 1; i < atomNum; i++) {
			r = getAtomicRadius(atomType[i]);
			rCube = r * r * r;
			xTotal += x[i] * rCube;
			yTotal += y[i] * rCube;
			zTotal += z[i] * rCube;
			rCubeTotal += rCube;
		}
		coreSite[1] = xTotal / rCubeTotal;
		coreSite[2] = yTotal / rCubeTotal;
		coreSite[3] = zTotal / rCubeTotal;
		return coreSite;
	}
	else {
		cerr << "Error: Invalid atom number or array sizes do not match atom number." << endl;
		coreSite[0] = -1;
		return coreSite;
	}
}

bool initArrayDect(int atomNum, std::vector<string> atomType, std::vector<double> x, std::vector<double> y, std::vector<double> z) {
	size_t xArraySize = x.size();
	size_t yArraySize = y.size();
	size_t zArraySize = z.size();
	size_t atomTypeSize = atomType.size();
	int atomNumReal = atomNum + 1;
	if (atomNumReal == xArraySize && atomNumReal == yArraySize && atomNumReal == zArraySize && atomNumReal == atomTypeSize) {
		return true;
	}
	else {
		return false;
	}
}