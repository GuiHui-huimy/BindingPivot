#include "../include/header.h"

void batchInstructions(string workPath) {
	cout << "-------------------------------------" << endl;
	cout << "You batch process work path is :" << workPath << endl;
	cout << "You need batchConfig.ini file in your work path." << endl;
	cout << "The more infomation ,Please check www.huimy.top ." << endl;
	cout << "-------------------------------------" << endl;
}

bool batchConfigInfo() {
	cout << "-------------------------------------" << endl;
	cout << "Batch Config File Instructions:" << endl;
	cout << "****************************************" << endl;
	cout << "haveLigandModel = " << getStringConfigValue("BATCH", "haveLigandModel", "") << endl;
	cout << "havePorcketModel = " << getStringConfigValue("BATCH", "havePorcketModel", "") << endl;
	cout << "****************************************" << endl;
	cout << "maxNumLigandModel = " << getStringConfigValue("BATCH", "maxNumLigandModel", "") << endl;
	cout << "maxNumPorcketModel = " << getStringConfigValue("BATCH", "maxNumPorcketModel", "") << endl;
	cout << "****************************************" << endl;
	cout << "colsedSiteCover = " << getStringConfigValue("BATCH", "colsedSiteCover", "") << endl;
	cout << "coverMaxDistance = " << getStringConfigValue("BATCH", "coverMaxDistance", "") << endl;
	cout << "****************************************" << endl;
	cout << "inputLigandModelName = " << getStringConfigValue("BATCH", "inputLigandModelName", "") << endl;
	cout << "inputPocketModelName = " << getStringConfigValue("BATCH", "inputPocketModelName", "") << endl;
	cout << "****************************************" << endl;
	cout << "ligandDirName = " << getStringConfigValue("BATCH", "ligandDirName", "") << endl;
	cout << "pocketDirName = " << getStringConfigValue("BATCH", "pocketDirName", "") << endl;
	cout << "****************************************" << endl;
	cout << "ligandProteinContentsFile = " << getStringConfigValue("BATCH", "ligandProteinContentsFile", "") << endl;
	cout << "pocketProteinContensFile = " << getStringConfigValue("BATCH", "pocketProteinContensFile", "") << endl;
	cout << "****************************************" << endl;
	cout << "-------------------------------------" << endl;
	cout << "Ok batch config file information ? (y/n) : ";
	char choice;
	cin >> choice;
	if (choice == 'y') {
		return true;
	}
	else {
		return false;
	}

}
ErrorHeader batchProcess(string workPath) {
	ErrorHeader ErrorInfo;
	int ErrorNum = 0;
	vector<string> ErrorList;

	batchInstructions(workPath);
	string batchConfigPath = workPath + "\\batchConfig.ini";
	ErrorInfo = readConfigData(batchConfigPath);
	if(ErrorInfo.isError){
		bool isCreatDefault = false;
		ErrorInfo = writeDefaultBatchConfig(workPath,isCreatDefault);
		if (isCreatDefault) {
			if (ErrorInfo.isError) {
				return ErrorInfo;
			}
		}
		else {
			cout << "please input your batchConfig.ini file in work path : " << endl;
			cin >> batchConfigPath;
		}
		ErrorInfo = readConfigData(batchConfigPath);
		if (ErrorInfo.isError) {
			ErrorInfo.errorMessage = "Error: Unable to open batch config file at :" + batchConfigPath;
			return ErrorInfo;
		}
	}
	if (!batchConfigInfo()) {
		ErrorInfo.isError = true;
		ErrorInfo.errorMessage = "Error: Batch config problem, process has been canceled by user.";
		return ErrorInfo;
	}


	bool haveLigandModel = getBoolConfigValue("BATCH", "haveLigandModel", false);
	bool havePorcketModel = getBoolConfigValue("BATCH", "havePorcketModel", false);

	if (!haveLigandModel && !havePorcketModel) {
		ErrorInfo.isError = true;
		ErrorInfo.errorMessage = "Error: In batch config file , both haveLigandModel and havePorcketModel are false. Please check your batchConfig.ini file.";
		return ErrorInfo;
	}
	if (haveLigandModel) {
		string ligandProteinListFilePath =workPath +"\\" + getStringConfigValue("BATCH", "ligandProteinContentsFile", "");
		vector<string> ligandProteinList;
		if (getProteinList(ligandProteinListFilePath, ligandProteinList)==SUCCESS) {
			
			if (ligandProteinList.size() == 0) {
				ErrorInfo.isError = true;
				ErrorInfo.errorMessage = "Error: In Ligand Protein ContentsFile File . The fi1e :" + ligandProteinListFilePath + " can't empty !";
				ErrorInfo.ErrorLogOut();
				ErrorNum++;
				ErrorList.push_back(ErrorInfo.errorMessage);
			} 
			else{
				cout << "Have found £º" << ligandProteinList.size() << "ligand analysis protein in list file" << endl;
				for (int i = 0; i < ligandProteinList.size(); i++) {
					trim(ligandProteinList[i]);
					string ligandModelpdbPath = workPath + "\\" + getStringConfigValue("BATCH", "ligandDirName", "") + "\\" + ligandProteinList[i] + "\\" + getStringConfigValue("BATCH", "inputLigandModelName", "model_");
					cout << "Start analysis ligand model protein :" << ligandProteinList[i] << endl;
					int maxLigandModelNum = getIntConfigValue("BATCH", "maxLigandModelNum", 3);
					bool isFirstError = true;
					for (int j = 1; j <= maxLigandModelNum; j++) {
						vector<double> returnCoreSite;
						vector <string> returnrecLigands;
						ErrorInfo = ligandBatchCore(ligandModelpdbPath + to_string(j) + ".pdb", ligandProteinList[i], getStringConfigValue("BATCH", "inputLigandModelName", "model_") + to_string(j), returnrecLigands, returnCoreSite);
						if (ErrorInfo.isError) {
							if (ErrorInfo.errorMessage == "Unable to open pdb file from :" + ligandModelpdbPath + to_string(j) + ".pdb") {
								if (isFirstError) {
									ErrorInfo.ErrorLogOut();
									ErrorNum++;
									ErrorList.push_back(ErrorInfo.errorMessage);
									string csvFilePath = workPath + "\\" + "CoreSiteBatchResults.csv";
									vector<string> initRows = { "Protein","CoreSite","Ligand or Pocket","Status","TestPdbLine" };
									vector<string> rows = { ligandProteinList[i],"","","Fail",""};

									if (writeCsvData(csvFilePath, initRows, rows) == ERROR) {
										ErrorInfo.isError = true;
										ErrorInfo.errorMessage = "Error: Unable to write csv file at :" + csvFilePath;
										return ErrorInfo;
									}
								}
								else {
									continue;
								}

							}
							else {
								ErrorInfo.ErrorLogOut();
								ErrorNum++;
								ErrorList.push_back(ErrorInfo.errorMessage);
								string csvFilePath = workPath + "\\" + "CoreSiteBatchResults.csv";
								vector<string> initRows = { "Protein","CoreSite","Ligand or Pocket","Status","TestPdbLine" };
								vector<string> rows = { ligandProteinList[i],"Ligand_Model","","Fail","" };

								if (writeCsvData(csvFilePath, initRows, rows) == ERROR) {
									ErrorInfo.isError = true;
									ErrorInfo.errorMessage = "Error: Unable to write csv file at :" + csvFilePath;
									return ErrorInfo;
								}
							}
						}
						else {
							isFirstError = false;
							int coreSiteIndex = 0;
							int ligandIndex = 0;
							vector<double> coreSiteSplit;
							for (int k = 0; k < returnCoreSite.size(); k++) {
								coreSiteSplit.push_back(returnCoreSite[k]);
								coreSiteIndex++;
								if (coreSiteIndex == 3) {
									if (coreSiteSplit.size() == 3) {
										string csvFilePath = workPath + "\\" + "CoreSiteBatchResults.csv";
										vector<string> initRows = { "Protein","CoreSite","Ligand or Pocket","Status","TestPdbLine" };
										string testPdbLine = createPdbLine(coreSiteSplit[0], coreSiteSplit[1], coreSiteSplit[2]);
										vector<string> rows = { ligandProteinList[i],formatDoubleToString(coreSiteSplit[0]) + "  " + formatDoubleToString(coreSiteSplit[1]) + "  " + formatDoubleToString(coreSiteSplit[2]),returnrecLigands[ligandIndex],"Success",testPdbLine};

										if (writeCsvData(csvFilePath, initRows, rows) == ERROR) {
											ErrorInfo.isError = true;
											ErrorInfo.errorMessage = "Error: Unable to write csv file at :" + csvFilePath;
											return ErrorInfo;
										}
										ligandIndex++;
										coreSiteIndex = 0;
										coreSiteSplit.clear();
									}

								}
							}
						}
					}

				}
			}
		}
		else {
			ErrorInfo.errorMessage = "Error: In Ligand Protein ContentsFile File . The fi1e :" + ligandProteinListFilePath + " can't empty !";
			ErrorInfo.ErrorLogOut();
			ErrorNum++;
			ErrorList.push_back(ErrorInfo.errorMessage);
		}
	}
	if (havePorcketModel) {
		string pocketProteinListFilePath = workPath + "\\" + getStringConfigValue("BATCH", "pocketProteinContensFile", "");
		vector<string> pocketProteinList;
		if (getProteinList(pocketProteinListFilePath,pocketProteinList) == SUCCESS) {
			if (pocketProteinList.size() == 0) {
				ErrorInfo.isError = true;
				ErrorInfo.errorMessage = "Error: In Ligand Protein ContentsFile File . The fi1e :" + pocketProteinListFilePath + " can't empty !";
				return ErrorInfo;
			}
			cout << "Have found £º" << pocketProteinList.size() << "pocket analysis protein in list file" << endl;
			for (int i = 0; i < pocketProteinList.size(); i++) {
				trim(pocketProteinList[i]);
				string pocketModelpdbPath = workPath + "\\" + getStringConfigValue("BATCH", "pocketDirName", "") + "\\" + pocketProteinList[i] + "\\" + getStringConfigValue("BATCH", "inputPocketModelName", "pocket_");
				cout << "Start analysis pocket model protein :" << pocketProteinList[i] << endl;
				int maxPocketModelNum = getIntConfigValue("BATCH", "maxNumPorcketModel", 3);
				bool isFirstError = true;
				int coreSiteIndex = 0;
				for (int j = 1; j <= maxPocketModelNum; j++) {
					vector<double> returnPocketCoreSite;
					int ligangNum = 1;
					string pocketProteinName = pocketProteinList[i] +"-" + getStringConfigValue("BATCH", "inputPocketModelName", "pocket_") + to_string(j);
					ErrorInfo = pocketBatchCore(pocketModelpdbPath + to_string(j) + ".pdb",pocketProteinName, getStringConfigValue("BATCH", "inputLigandModelName", "pocket_") + to_string(j), returnPocketCoreSite);
					if (ErrorInfo.isError) {
						if (ErrorInfo.errorMessage == "Unable to open pdb file from : " + pocketModelpdbPath + to_string(j) + ".pdb") {
							if (isFirstError) {
								ErrorInfo.ErrorLogOut();
								ErrorNum++;
								ErrorList.push_back(ErrorInfo.errorMessage);
								string csvFilePath = workPath + "\\" + "CoreSiteBatchResults.csv";
								vector<string> initRows = { "Protein","CoreSite","Ligand or Pocket","Status","TestPdbLine" };
								vector<string> rows = { pocketProteinList[i],"","","Fail","" };

								if (writeCsvData(csvFilePath, initRows, rows) == ERROR) {
									ErrorInfo.isError = true;
									ErrorInfo.errorMessage = "Error: Unable to write csv file at :" + csvFilePath;
									return ErrorInfo;
								}
							}
							else {
								continue;
							}
						}
						else {
							ErrorInfo.ErrorLogOut();
							ErrorNum++;
							ErrorList.push_back(ErrorInfo.errorMessage);
							string csvFilePath = workPath + "\\" + "CoreSiteBatchResults.csv";
							vector<string> initRows = { "Protein","CoreSite","Ligand or Pocket","Status","TestPdbLine" };
							vector<string> rows = { pocketProteinList[i],"Pocket_Model","","Fail","" };

							if (writeCsvData(csvFilePath, initRows, rows) == ERROR) {
								ErrorInfo.isError = true;
								ErrorInfo.errorMessage = "Error: Unable to write csv file at :" + csvFilePath;
								return ErrorInfo;
							}
						
						
						
						}
					}
					else {
						isFirstError = false;
						int coreSiteIndex = 0;
						vector<double> coreSiteSplitPocket;
						for (int k = 0; k < returnPocketCoreSite.size(); k++) {
							coreSiteSplitPocket.push_back(returnPocketCoreSite[k]);
							coreSiteIndex++;

							if (coreSiteIndex == 3) {
								if (coreSiteSplitPocket.size() == 3) {
									string csvFilePath = workPath + "\\" + "CoreSiteBatchResults.csv";
									vector<string> initRows = { "Protein","CoreSite","Ligand or Pocket","Status","TestPdbLine" };
									string testPdbLine = createPdbLine(coreSiteSplitPocket[0], coreSiteSplitPocket[1], coreSiteSplitPocket[2]);
									//vector<string> rows = {pocketProteinList[i],to_string(coreSiteSplitPocket[0]) + "  " + to_string(coreSiteSplitPocket[1]) + "  " + to_string(coreSiteSplitPocket[2]),testPdbLine};
									vector<string> rows = { pocketProteinList[i],formatDoubleToString(coreSiteSplitPocket[0]) + "  " + formatDoubleToString(coreSiteSplitPocket[1]) + "  " + formatDoubleToString(coreSiteSplitPocket[2]),getStringConfigValue("BATCH", "inputPocketModelName", "pocket_") + to_string(j),"Success",testPdbLine };
									if (writeCsvData(csvFilePath, initRows, rows) == ERROR) {
										ErrorInfo.isError = true;
										ErrorInfo.errorMessage = "Error: Unable to write csv file at :" + csvFilePath;
										return ErrorInfo;
									}
									coreSiteIndex = 0;
									coreSiteSplitPocket.clear();
								}

							}
						}
					}
				}
			}
		}
	}

	cout << "============================================" << endl;
	cout << "Batch process completed. Total proteins processed: " << (haveLigandModel ? getIntConfigValue("BATCH", "maxLigandModelNum", 3) * getIntConfigValue("BATCH", "maxNumPorcketModel", 3) : 0) + (havePorcketModel ? getIntConfigValue("BATCH", "maxNumPorcketModel", 3) * getIntConfigValue("BATCH", "maxNumPorcketModel", 3) : 0) << endl;
	cout << "Batch process completed with Non-fatal Error £º" << ErrorNum << " errors." << endl;
	if (ErrorNum > 0) {
		cout << "Error List:" << endl;
		for (int i = 0; i < ErrorList.size(); i++) {
			cout << i + 1 << ". " << ErrorList[i] << endl;
		}
	}
	else {
		cout << "No errors encountered during batch process." << endl;
	}
	ErrorInfo.isError = false;
	return ErrorInfo;
}

