#include "../include/header.h"

ErrorHeader IndexPage(string& workPath) {
	cout << "Welcome to AOTM Core Site Project!" << endl;
	cout << "This project is used to calculate  core site based on pdb file data." << endl;
	cout << "Now this version is :" << getStringConfigValue("GENERAL", "version", "1.00")<<endl;
	ErrorHeader ErrorInfo = setWorkPath(workPath);
	return ErrorInfo;
}
ErrorHeader setWorkPath(string& workPath) {
	ErrorHeader ErrorInfo;
	string defaultWorkPath = getStringConfigValue("FILE", "defaultWorkPath", "");
Reset:
	cout << "-------------------------------------" << endl;
	cout << "Please set workpath to AOTM Core Site project." << endl;
	if (!defaultWorkPath.empty()) {
		cout << "If you input empty will use default work path. " << endl;
		cout <<"Default work path:" << defaultWorkPath << endl;
	}else {
		cout << "You can set default work path in config file." << endl;
	}
	cout << "WorkPath dir : ";
	cin >> workPath;
	if (!workPath.empty()) {
		ErrorInfo.isError = false;
		return ErrorInfo;
	}
	else {
		if(!defaultWorkPath.empty()){
			workPath = defaultWorkPath;
			ErrorInfo.isError = false;
			return ErrorInfo;
		}
		else {
			cout << "None work path !\n";
			goto Reset;
		}
	}

}
ErrorHeader chosePdbAnalysis(string& workPath ,string& fullPdbPath) {
	ErrorHeader ErrorInfo;
	string pdbFileName = "";
	fullPdbPath = "";
	Reinput_pdbname:
	cout << "Please input you .pdb file name. Just like : model_1.pdb " << endl;
	cout << "File name : ";
	cin >> pdbFileName;
	if (!pdbFileName.empty()) {
		fullPdbPath = workPath + "\\" + pdbFileName;
		ErrorInfo.isError = false;
		return ErrorInfo;
	}
	else {
		cout << "None pdb file name !\n";
		goto Reinput_pdbname;
	}

}
ErrorHeader choseList(string& workPath) {
	string fullPdbPath = "";
	int choice = 0;
	ErrorHeader ErrorInfo;
	Reinput_choice:
	cout << "-------------------------------------" << endl;
	cout << "Reading PDB file workpath: " << workPath << endl;
	cout << "Please chose Action you want to do.\n" << endl;
	cout << "**********************************" << endl;
	cout << "1.Reset workPath" << endl;
	cout << "2.Analysis ligand core site from complex pdb file" << endl;
	cout << "3.Analysis residue from protein pdb file" << endl;
	cout << "4.Analysis pocket residue from pocket pdb file" << endl;
	cout << "5.Batch analysis core site" << endl;
	cout << "**********************************\n" << endl;
	cout << "Enter choice:";
	cin >> choice;
	if (choice == 1 || choice == 2||choice == 3 ||choice == 4||choice == 5) {
		switch (choice) {
		case 1:
			ErrorInfo = setWorkPath(workPath);
			if (ErrorInfo.isError) { return ErrorInfo; }
			break;
		case 2:
			cout << "-------------------------------------" << endl;
			cout << "Chose complex PDB file analysis Ligand core site from workPath: " << workPath << endl;
			ErrorInfo = chosePdbAnalysis(workPath,fullPdbPath);
			if (ErrorInfo.isError) { return ErrorInfo; }
			ErrorInfo = analysisLigandCoreSite(fullPdbPath);
			if (ErrorInfo.isError) { return ErrorInfo; }
			break;
		case 3:
			cout << "-------------------------------------" << endl;
			cout << "Chose protein PDB file analysis Res core site from workPath: " << workPath << endl;
			ErrorInfo = chosePdbAnalysis(workPath,fullPdbPath);
			if (ErrorInfo.isError) { return ErrorInfo; }
			ErrorInfo = analysisResCoreSite(fullPdbPath);
			if (ErrorInfo.isError) { return ErrorInfo; }
			break;
		case 4:
			cout << "-------------------------------------" << endl;
			cout << "Chose pocket PDB file analysis Pocket core site from workPath: " << workPath << endl;
			ErrorInfo = chosePdbAnalysis(workPath, fullPdbPath);
			if(ErrorInfo.isError){return ErrorInfo;}
			ErrorInfo =	analysisPocketCoreSite(fullPdbPath);
			if (ErrorInfo.isError) { return ErrorInfo; }
			break;
		case 5:
			cout << "-------------------------------------" << endl;
			cout << "Chose batch analysis core site from workPath: " << workPath << endl;
			ErrorInfo = batchProcess(workPath);
			if (ErrorInfo.isError) { return ErrorInfo; }
			break;

		}
		return ErrorInfo;
	}
	else {
		cout << "Invalid choice ! Please chose again.\n";
		goto Reinput_choice;
	}
}