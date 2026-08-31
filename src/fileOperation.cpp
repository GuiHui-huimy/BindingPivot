#include "../include/header.h"

void trim(std::string& str) {
    str.erase(0, str.find_first_not_of(" \t\n\r\f\v"));
    str.erase(str.find_last_not_of(" \t\n\r\f\v") + 1);
}

int getPdbData(string fileDir, vector<string>& pdbLines) {
    std::ifstream pdbFile(fileDir);
    if (pdbFile.is_open()) {
        string line;
        while (std::getline(pdbFile, line))
        {
            pdbLines.push_back(line);
        }
        pdbFile.close();
        return SUCCESS;
    }
    else {
        return ERROR;

    }

}

int getConfigData(string fileDir, vector<string>& configLines) {
    std::ifstream configFile(fileDir);
    if (configFile.is_open()) {
        string line;
        while (std::getline(configFile, line))
        {
            configLines.push_back(line);
        }
        configFile.close();
        return SUCCESS;
    }
    else {
        return ERROR;

    }

}

int  getProteinList(string fileDir,vector<string>&  proteinList ) {
    std::ifstream listFile(fileDir);
    if (listFile.is_open()) {
        string line;
        while (std::getline(listFile, line)) {
            proteinList.push_back(line);
        }
		listFile.close();
		return SUCCESS;
    }
    else {
		return ERROR;
    }
}

int writeCsvData(string fileDir,vector<string> initRows,vector<string> rows) {

    std::ifstream iCsvFile(fileDir);
    if (iCsvFile.good()) {
		iCsvFile.close();
		std::ofstream oCsvFile(fileDir, std::ios::app);
        if (oCsvFile.is_open()) {
            for (int i = 0; i < rows.size();i++) {
				oCsvFile << rows[i];
                if (i != rows.size() - 1) {
					oCsvFile << ",";
                }
            }
			oCsvFile << "\n";
			oCsvFile.close();
			return SUCCESS;
        }
        else
        {
			return ERROR;
        }
    }
    else {
		iCsvFile.close();
        std::ofstream oCsvFile(fileDir);
        if (oCsvFile.is_open()) {
            for (int i = 0; i < initRows.size(); i++) {
                oCsvFile << initRows[i];
                if (i != initRows.size() - 1) {
                    oCsvFile << ",";
                }
            }
            oCsvFile << "\n";
            for (int i = 0; i < rows.size(); i++) {
                oCsvFile << rows[i];
                if (i != rows.size() - 1) {
                    oCsvFile << ",";
                }
            }
            oCsvFile << "\n";
			oCsvFile.close();
			return SUCCESS;
        }
        else
        {
			return ERROR;
        }




    }

}
ErrorHeader writeDefaultBatchConfig(string workPath ,bool& creatDefault) {
	cout << "No batchConfig.ini file found in work path : " << workPath << endl;
	cout << "Would you like to create a default batchConfig.ini file in this path ? (Y/N) : " << endl;
	char userInput;
	cin >> userInput;
	ErrorHeader ErrorInfo;
    if (userInput == 'Y' || userInput == 'y') {
        creatDefault = true;
        std::ofstream configFile(workPath + "\\" + "batchConfig.ini");
        if (configFile.is_open()) {
            configFile << "[BATCH]\n";
            configFile << "\n";
            configFile << "#批量是否有配体模型\n";
            configFile << "haveLigandModel = true\n";
            configFile << "#批量是否有口袋模型\n";
            configFile << "havePorcketModel = true\n";
            configFile << "\n";
            configFile << "#最大配体模型数量\n";
            configFile << "maxNumLigandModel = 3\n";
            configFile << "#最大口袋模型数量\n";
            configFile << "maxNumPorcketModel = 3\n";
            configFile << "\n";
            configFile << "#空间接近位点覆盖\n";
            configFile << "colsedSiteCover = true\n";
            configFile << "#覆盖检测阈值：一般是对接口袋立方体边长\n";
            configFile << "coverMaxDistance = 20\n";
            configFile << "\n";
            configFile << "#配体文件序列前缀名称：例如\"model-1\"的前缀名称为\"model-\"\n";
            configFile << "inputLigandModelName = model_\n";
            configFile << "#口袋文件序列前缀名称：例如\"pocket-1\"的前缀名称为\"pocket-\"\n";
            configFile << "inputPocketModelName = pocket_\n";
            configFile << "\n";
            configFile << "#配体批量文件夹名称\n";
            configFile << "ligandDirName = ligand\n";
            configFile << "#口袋批量文件夹名称\n";
            configFile << "pocketDirName = pocket\n";
            configFile << "\n";
            configFile << "#存放需要由配体分析蛋白质名称汇总的文件\n";
            configFile << "ligandProteinContentsFile = ligandProtein.txt\n";
            configFile << "#存放需要由口袋分析蛋白质名称汇总文件\n";
            configFile << "pocketProteinContensFile = pocketProtein.txt\n";
            configFile.close();
			ErrorInfo.isError = false;
			return ErrorInfo;
        }
        else {
			ErrorInfo.isError = true;
			ErrorInfo.errorMessage = "Error: Unable to create default batchConfig.ini file in path :" + workPath;
			return ErrorInfo;
 
        }
    }
    else {
		creatDefault = false;
		ErrorInfo.isError = false;
		return ErrorInfo;
    }
}
