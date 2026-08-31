#include "../include/header.h"
ErrorHeader analysisLigandCoreSite(string fileDir) {
	//initialize variables
	ErrorHeader ErrorInfo;
	int i;
	std::vector<double> x(0);
	std::vector<double> y(0);
	std::vector<double> z(0);
	std::vector<string> atomType;
	int atomNum = 1;
	std::vector<double> coreSite(4);
	string inputLines;
	vector<string> Lines;
	bool isHaveThisMol = false;
	vector <string> recLigands;
	vector <int> recLigandseries;
	size_t lineNum;

	if (initGetInfoFromLines(lineNum,fileDir,Lines,x,y,z,atomType)==SUCCESS) {
		
		recogLigand(lineNum, Lines, recLigands,recLigandseries);
		if (!recLigands.empty()&&!recLigandseries.empty()) {
			if (recLigands.size() == recLigandseries.size()) {
				bool showTableBar = true;
				showTableBar = getBoolConfigValue("DISPLAY", "showTableBar", true);
				if (showTableBar) {
					cout << "©°©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©Ð©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©´\n";
					cout << "©¦ Ligand Name      ©¦ Serial Number   ©¦\n";
					cout << "©À©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©à©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©È\n";
					for (size_t i = 0; i < recLigands.size(); i++) {
						cout << "©¦ " << setw(16) << left << recLigands[i]
							<< " ©¦ " << setw(15) << recLigandseries[i] << " ©¦\n";
					}
					cout << "©¸©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©Ø©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¼\n";
				}
				else {
					cout << "**********************************" << endl;
					cout << "Recognized Ligand Molecules in PDB file:" << endl;
					for (size_t i = 0; i < recLigands.size(); i++)
					{
						cout << "Ligand Name   Serial Number" << "\n" << endl;
						cout << recLigands[i] << "           " << recLigandseries[i] << "\n" << endl;
						cout << "**********************************" << endl;
					}
				}
			
			}
			else
			{
				ErrorInfo.isError = true;
				ErrorInfo.errorMessage = "Ligand Molecules recognized number and serial number size not match.";
				return ErrorInfo;
			}
		
		}
		else {
			ErrorInfo.isError = true;
			ErrorInfo.errorMessage = "No Ligand Molecules recognized in PDB file.";
			return ErrorInfo;
		}
	}
	else
	{
		ErrorInfo.isError = true;
		ErrorInfo.errorMessage = "Unable to open pdb file from :" + fileDir;
		return ErrorInfo;
	}
	
	cout << "Please input ligand molecule name:";
	string ligandName;
	cin >> ligandName;
	if (ligandName.empty()) {
		ErrorInfo.isError = true;
		ErrorInfo.errorMessage = "Ligand name can't empty";
		return ErrorInfo;
	}
	cout << "Please input ligand molecule  serial Number:";
	int serialNum = 0;
	cin >> serialNum;
	if (serialNum == 0) {
		ErrorInfo.isError = true;
		ErrorInfo.errorMessage = "Ligand name can't empty.";
		return ErrorInfo;
	}

	bool showLinesInfo = true;
	bool showWarnings = true;
	bool showSpecialLine = true;
	showLinesInfo = getBoolConfigValue("DISPLAY", "showLinesInfo", true);
	showWarnings = getBoolConfigValue("DISPLAY", "showWarnings", true);
	showSpecialLine = getBoolConfigValue("DISPLAY", "showSpecialLine", true);

	// Loop to get atom information
	for (i = 1; i <= lineNum; i++) {
		cout << "Read line - " << i << "Info " << endl;
		AtomInfo atomInfos;
		if (specialLine(Lines[i - 1], showSpecialLine)) { break; }
		if (extractSingleLine(Lines[i-1], atomInfos)) {
			if (isNeedMol(ligandName, serialNum, atomInfos)==true) {
				if (showLinesInfo) {
					cout << "Ligand Atom Type: " << atomInfos.type << ", Coordinates: (" << atomInfos.x << ", " << atomInfos.y << ", " << atomInfos.z << ")" << endl;
				}
				x[atomNum] = atomInfos.x;
				y[atomNum] = atomInfos.y;
				z[atomNum] = atomInfos.z;
				atomType[atomNum] = atomInfos.type;
				isHaveThisMol = true;
				atomNum++;
			}
		}
		else {
			if(showWarnings){
				ErrorInfo.isError = true;
				ErrorInfo.errorMessage = "Invalid atom information format. Error in line "+ i ;
				return ErrorInfo;
			}
		}

	}
	if (isHaveThisMol == false) {
		ErrorInfo.isError = true;
		ErrorInfo.errorMessage = "Can't find this ligand Molecule:" + ligandName;
		return ErrorInfo;
	}



	cout << "-------------------------------------" << endl;
	// Calculate core site
	atomType.resize(atomNum + 1);
	x.resize(atomNum + 1);
	y.resize(atomNum + 1);
	z.resize(atomNum + 1);
	coreSite = rWeightedSiteCalc(atomNum, atomType, x, y, z);
	cout << "The ligand Molecule :" << ligandName << endl;
	cout << "Core Site Coordinates: (" << coreSite[0] << ", " << coreSite[1] << ", " << coreSite[2] << "," << coreSite[3] << ")" << endl;
	cout << fixed << setprecision(3) << "HETATM   33  COR MIX     1       " << coreSite[1] << "   " << coreSite[2] << "  " << coreSite[3] << " - TestPdbLine" << endl;
	cout << "Welcome to the AOTM Core Site Next!" << endl;
	ErrorInfo.isError = false;
	return ErrorInfo;
}

ErrorHeader analysisResCoreSite(string fileDir) {
	//initialize variables
	ErrorHeader ErrorInfo;
	int i;
	std::vector<double> x(0);
	std::vector<double> y(0);
	std::vector<double> z(0);
	std::vector<string> atomType;
	int atomNum = 1; 
	std::vector<double> coreSite(4);
	string inputLines;
	vector<string> Lines;
	bool isHaveThisMol = false;
	rsize_t lineNum;

	if (initGetInfoFromLines(lineNum, fileDir, Lines, x, y, z, atomType) == ERROR) {
		ErrorInfo.isError = true;
		ErrorInfo.errorMessage = "Unable to open pdb file from :" + fileDir;
		return ErrorInfo;
	}
	cout << "Please input Res molecule name:";
	string resName;
	cin >> resName;
	if (resName.empty()) {
		ErrorInfo.isError = true;
		ErrorInfo.errorMessage = "Res name can't empty.";
		return ErrorInfo;
	}
	cout << "Please input res molecule  serial Number:";
	int serialNum = 0;
	cin >> serialNum;
	if (serialNum == 0) {
		ErrorInfo.isError = true;
		ErrorInfo.errorMessage = "Res serial can't empty.";
		return ErrorInfo;
	}
	bool showLinesInfo = true;
	bool showWarnings = true;
	bool showSpecialLine = true;
	showLinesInfo = getBoolConfigValue("DISPLAY", "showLinesInfo", true);
	showWarnings = getBoolConfigValue("DISPLAY", "showWarnings", true);
	showSpecialLine = getBoolConfigValue("DISPLAY", "showSpecialLine", true);
	// Loop to get atom information
	for (i = 1; i <= lineNum; i++) {
		cout << "Read line - " << i << "Info " << endl;
		AtomInfo atomInfos;
		if (specialLine(Lines[i - 1], showSpecialLine)) { break; }
		if (extractSingleLine(Lines[i - 1], atomInfos)) {
			if (isNeedMol(resName, serialNum, atomInfos) == true) {
				if (showLinesInfo) {
					cout << "Res Atom Type: " << atomInfos.type << ", Coordinates: (" << atomInfos.x << ", " << atomInfos.y << ", " << atomInfos.z << ")" << endl;
				}
				x[atomNum] = atomInfos.x;
				y[atomNum] = atomInfos.y;
				z[atomNum] = atomInfos.z;
				atomType[atomNum] = atomInfos.type;
				isHaveThisMol = true;
				atomNum++;
			}
		}
		else {
			if(showWarnings){
				ErrorInfo.isError = true;
				ErrorInfo.errorMessage = "Invalid atom information format. Error in line " + i;
				return ErrorInfo;

			}
		}

	}
	if (isHaveThisMol == false) {
		ErrorInfo.isError = true;
		ErrorInfo.errorMessage = "Can't find this res Molecule: " + resName;
		return ErrorInfo;
	}
	cout << "-------------------------------------" << endl;
	// Calculate core site
	atomType.resize(atomNum + 1);
	x.resize(atomNum + 1);
	y.resize(atomNum + 1);
	z.resize(atomNum + 1);
	coreSite = rWeightedSiteCalc(atomNum, atomType, x, y, z);
	cout << "The res Molecule :" << resName << endl;
	cout << "Core Site Coordinates: (" << coreSite[0] << ", " << coreSite[1] << ", " << coreSite[2] << "," << coreSite[3] << ")" << endl;
	cout << fixed << setprecision(3) << "HETATM   33  COR MIX     1       " << coreSite[1] << "   " << coreSite[2] << "  " << coreSite[3] << " - TestPdbLine" << endl;
	cout << "Welcome to the AOTM Core Site Next!" << endl;
	ErrorInfo.isError = false;
	return ErrorInfo;
}
ErrorHeader analysisPocketCoreSite(string fileDir) {
	//initialize variables
	ErrorHeader ErrorInfo;
	int i;
	std::vector<double> x(0);
	std::vector<double> y(0);
	std::vector<double> z(0);
	std::vector<string> atomType;
	int atomNum = 1;
	std::vector<double> coreSite(4);
	string inputLines;
	vector<string> Lines;
	rsize_t lineNum;

	if (initGetInfoFromLines(lineNum, fileDir, Lines, x, y, z, atomType) == ERROR) {
		ErrorInfo.isError = true;
		ErrorInfo.errorMessage = "Unable to open pdb file from :" + fileDir;
		return ErrorInfo;
	}
	bool showLinesInfo = true;
	bool showWarnings = true;
	bool showSpecialLine = true;
	showLinesInfo = getBoolConfigValue("DISPLAY", "showLinesInfo", true);
	showWarnings = getBoolConfigValue("DISPLAY", "showWarnings", true);
	showSpecialLine = getBoolConfigValue("DISPLAY", "showSpecialLine", true);

	// Loop to get atom information
	for (i = 1; i <= lineNum; i++) {
		cout << "Read line - " << i << "Info " << endl;
		AtomInfo atomInfos;
		if (specialLine(Lines[i - 1], showSpecialLine)) { break; }
		if (extractSingleLine(Lines[i - 1], atomInfos)) {
			if (showLinesInfo) {
				cout << "Pocket Atom Type: " << atomInfos.type << ", Coordinates: (" << atomInfos.x << ", " << atomInfos.y << ", " << atomInfos.z << ")" << endl;
			}
			x[atomNum] = atomInfos.x;
			y[atomNum] = atomInfos.y;
			z[atomNum] = atomInfos.z;
			atomType[atomNum] = atomInfos.type;
			atomNum++;
		}
		else {
			if (showWarnings) {
				ErrorInfo.isError = true;
				ErrorInfo.errorMessage = "Invalid atom information format. Error in line " + i;
				return ErrorInfo;
			}
		}

	}
	cout << "-------------------------------------" << endl;
	// Calculate core site
	atomType.resize(atomNum + 1);
	x.resize(atomNum + 1);
	y.resize(atomNum + 1);
	z.resize(atomNum + 1);
	coreSite = rWeightedSiteCalc(atomNum, atomType, x, y, z);
	cout << "The Pocket :" << fileDir << endl;
	cout << "Core Site Coordinates: (" << coreSite[0] << ", " << coreSite[1] << ", " << coreSite[2] << "," << coreSite[3] << ")" << endl;
	cout << fixed << setprecision(3) << "HETATM   33  COR MIX     1       " << coreSite[1] << "   " << coreSite[2] << "  " << coreSite[3] << " - TestPdbLine" << endl;
	cout << "Welcome to the AOTM Core Site Next!" << endl;
	ErrorInfo.isError =false;
	return ErrorInfo;
}

ErrorHeader ligandBatchCore(string fileDir, string ligandProteinName, string modelName, vector <string>& returnrecLigands, vector<double>& returnCoreSite) {
	//initialize variables
	ErrorHeader ErrorInfo;

	std::vector<double> x(0);
	std::vector<double> y(0);
	std::vector<double> z(0);
	std::vector<string> atomType;
	int atomNum = 1;
	std::vector<double> coreSite(4);
	string inputLines;
	vector<string> Lines;
	vector <string> recLigands;
	vector <int> recLigandseries;
	size_t lineNum;

	if (initGetInfoFromLines(lineNum, fileDir, Lines, x, y, z, atomType) == SUCCESS) {

		recogLigand(lineNum, Lines, recLigands, recLigandseries);
		if (!recLigands.empty() && !recLigandseries.empty()) {
			if (recLigands.size() == recLigandseries.size()) {
				returnrecLigands = recLigands;
				bool showTableBar = true;
				showTableBar = getBoolConfigValue("DISPLAY", "showTableBar", true);
				if (showTableBar) {
					cout << "©°©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©Ð©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©´\n";
					cout << "©¦ Ligand Name      ©¦ Serial Number   ©¦\n";
					cout << "©À©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©à©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©È\n";
					for (size_t i = 0; i < recLigands.size(); i++) {
						cout << "©¦ " << setw(16) << left << recLigands[i]
							<< " ©¦ " << setw(15) << recLigandseries[i] << " ©¦\n";
					}
					cout << "©¸©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©Ø©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¼\n";
				}
				else {
					cout << "**********************************" << endl;
					cout << "Recognized Ligand Molecules in PDB file:" << endl;
					for (size_t i = 0; i < recLigands.size(); i++)
					{
						cout << "Ligand Name   Serial Number" << "\n" << endl;
						cout << recLigands[i] << "           " << recLigandseries[i] << "\n" << endl;
						cout << "**********************************" << endl;
					}
				}

			}
			else
			{
				ErrorInfo.isError = true;
				ErrorInfo.errorMessage = "Ligand Molecules recognized number and serial number size not match. In protein :" + ligandProteinName + " In ligand model :" + modelName;
				return ErrorInfo;
			}

		}
		else {
			ErrorInfo.isError = true;
			ErrorInfo.errorMessage = "No Ligand Molecules recognized in PDB file. In file : " + fileDir;
			return ErrorInfo;
		}
	}
	else
	{
		ErrorInfo.isError = true;
		ErrorInfo.errorMessage = "Unable to open pdb file from :" + fileDir;
		return ErrorInfo;
	}
	// Loop to get atom information
	for (int j = 0; j < recLigands.size(); j++) {
		atomNum = 1;
		for (int i = 1; i <= lineNum; i++) {
			//cout << "Read line - " << i << "Info " << endl;
			AtomInfo atomInfos;
			if (specialLine(Lines[i - 1], false)) { continue; }
			if (extractSingleLine(Lines[i - 1], atomInfos)) {
				if (isNeedMol(recLigands[j], recLigandseries[j], atomInfos) == true) {
					x[atomNum] = atomInfos.x;
					y[atomNum] = atomInfos.y;
					z[atomNum] = atomInfos.z;
					atomType[atomNum] = atomInfos.type;
					atomNum++;
				}
			}
			else {
				ErrorInfo.isError = true;
				ErrorInfo.errorMessage = ligandProteinName + ": Invalid atom information format. Error in line " + to_string(i);
				return ErrorInfo;

			}
		}
		cout << "-------------------------------------" << endl;
		// Calculate core site
		atomType.resize(atomNum + 1);
		x.resize(atomNum + 1);
		y.resize(atomNum + 1);
		z.resize(atomNum + 1);
		coreSite = rWeightedSiteCalc(atomNum, atomType, x, y, z);
		cout << "The ligand :" << recLigands[j] << endl;
		cout << "Core Site Coordinates: (" << coreSite[0] << ", " << coreSite[1] << ", " << coreSite[2] << "," << coreSite[3] << ")" << endl;
		cout << fixed << setprecision(3) << "HETATM   33  COR MIX     1       " << coreSite[1] << "   " << coreSite[2] << "  " << coreSite[3] << " - TestPdbLine" << endl;
		returnCoreSite.push_back(coreSite[1]);
		returnCoreSite.push_back(coreSite[2]);
		returnCoreSite.push_back(coreSite[3]);
	}
	ErrorInfo.isError = false;
	return ErrorInfo;
}

ErrorHeader pocketBatchCore(string fileDir, string pocketProteinName, string modelName, vector<double>& returnCoreSite) {
	//initialize variables
	ErrorHeader ErrorInfo;
	int i;
	std::vector<double> x(0);
	std::vector<double> y(0);
	std::vector<double> z(0);
	std::vector<string> atomType;
	int atomNum = 1;
	std::vector<double> coreSite(4);
	string inputLines;
	vector<string> Lines;
	vector <string> recLigands;

	size_t lineNum;

	if (initGetInfoFromLines(lineNum, fileDir, Lines, x, y, z, atomType) != SUCCESS) {
		ErrorInfo.isError = true;
		ErrorInfo.errorMessage = "Unable to open pdb file from :" + fileDir;
		return ErrorInfo;
	}
	// Loop to get atom information

	atomNum = 1;
	for (i = 1; i <= lineNum; i++) {
		//cout << "Read line - " << i << "Info " << endl;
		AtomInfo atomInfos;
		if (specialLine(Lines[i - 1], false)) { continue; }
		if (extractSingleLine(Lines[i - 1], atomInfos)) {
			x[atomNum] = atomInfos.x;
			y[atomNum] = atomInfos.y;
			z[atomNum] = atomInfos.z;
			atomType[atomNum] = atomInfos.type;
			atomNum++;
		}
		else {
			ErrorInfo.isError = true;
			ErrorInfo.errorMessage = "Invalid atom information format. Error in line " + i;
			return ErrorInfo;
		}
	}
	cout << "-------------------------------------" << endl;
	// Calculate core site
	atomType.resize(atomNum + 1);
	x.resize(atomNum + 1);
	y.resize(atomNum + 1);
	z.resize(atomNum + 1);
	coreSite = rWeightedSiteCalc(atomNum, atomType, x, y, z);
	cout << "The pocket :" << pocketProteinName << endl;
	cout << "Core Site Coordinates: (" << coreSite[0] << ", " << coreSite[1] << ", " << coreSite[2] << "," << coreSite[3] << ")" << endl;
	cout << fixed << setprecision(3) << "HETATM   33  COR MIX     1       " << coreSite[1] << "   " << coreSite[2] << "  " << coreSite[3] << " - TestPdbLine" << endl;
	returnCoreSite.push_back(coreSite[1]);
	returnCoreSite.push_back(coreSite[2]);
	returnCoreSite.push_back(coreSite[3]);

	ErrorInfo.isError = false;
	return ErrorInfo;
}




ErrorHeader recogLigand(int lineNum, vector<string> Lines, vector <string>& recLigands ,vector<int>& recLigandseries) {
	ErrorHeader ErrorInfo;
	int i = 1;
	
	for (i = 1; i <= lineNum; i++) {
		AtomInfo atomInfos;
		if (specialLine(Lines[i - 1], false)) { break; }
		if (extractSingleLine(Lines[i - 1], atomInfos)) {
			if (atomInfos.recordName == "HETATM") {
				if (std::find(recLigands.begin(), recLigands.end(), atomInfos.moleculeName) != recLigands.end()) {
					if (std::find(recLigands.begin(), recLigands.end(), atomInfos.moleculeName) != recLigands.end()) {
						break;
					}
					else {
						recLigands.push_back(atomInfos.moleculeName);
						recLigandseries.push_back(atomInfos.moleculeSerial);
					}
					
				}
				else {
					recLigands.push_back(atomInfos.moleculeName);
					recLigandseries.push_back(atomInfos.moleculeSerial);

				}
			}
		}
		else {
			ErrorInfo.isError = true;
			ErrorInfo.errorMessage = "Invalid atom information format. Error in line " + i;
			return ErrorInfo;
		}

	}
	ErrorInfo.isError = false;
	return ErrorInfo;
}
bool isNeedMol(string MolName, int MolSer, AtomInfo infos) {
	if (infos.moleculeSerial == MolSer && infos.moleculeName == MolName) {
		return true;
	}
	else {
		return false;
	}
}