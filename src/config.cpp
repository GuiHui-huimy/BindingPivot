#include "../include/header.h"
#include <windows.h>
map<string, map<string, string>> configData;
string getExePath() {
    char buffer[MAX_PATH];
    GetModuleFileNameA(NULL, buffer, MAX_PATH);
    std::string::size_type pos = std::string(buffer).find_last_of("\\/");
    return std::string(buffer).substr(0, pos);
}


ErrorHeader loadConfig() {
    ErrorHeader ErrorInfo;
    string exeDir = getExePath();
    string generalConfigPath = exeDir + "\\config\\general.ini";
    string analysisConfigPath = exeDir + "\\config\\analysis.ini";
    ErrorInfo = readConfigData(generalConfigPath);
    if (ErrorInfo.isError) { return ErrorInfo; }
    ErrorInfo = readConfigData(analysisConfigPath);
    if (ErrorInfo.isError) { return ErrorInfo; }
	return ErrorInfo;
}


ErrorHeader readConfigData(string ConfigPath) {
    ErrorHeader ErrorInfo;
    vector<string> configLines;
    if (!getConfigData(ConfigPath, configLines)) {
        ErrorInfo.isError = true;
        ErrorInfo.errorMessage = "Error: Unable to open config file at" + ConfigPath;
        return ErrorInfo;
    }
	string currentSection = "";
    for (int i = 0; i < configLines.size(); i++) {
        if (isCommendLine(configLines[i])) {
            continue;
        }
        if (isSectionHeader(configLines[i], currentSection)) {
			configData[currentSection] = map<string, string>();
			currentSection = currentSection;
            continue;
        }
        string  key, value;
        if (!configLines[i].empty()) {
            if (pharseKeyValue(configLines[i], key, value)) {
                configData[currentSection][key] = value;
            }
        }
        else {
            continue;
        }
    }
    ErrorInfo.isError = false;
    return ErrorInfo;
}

bool isSectionHeader(string configLine, string& sectionName) {
    if (!configLine.empty()) {
        string trimed = configLine;
        trim(trimed);
        if (configLine[0] == '[' && configLine.back() == ']') {
            sectionName = trimed.substr(1,trimed.length() - 2);
            return true;
        }
        else {
            return false;
        }
    }
    else {
        return false;
    }
}

bool isCommendLine(string configLine) {
        string trimed = configLine;
        trim(trimed);
        if (trimed[0] == '#' || trimed[0] == ';' || trimed[0]=='//') {
            return true;
        }
        else {
            return false;
        }
}

bool pharseKeyValue(string configLine, string& key, string& value) {
    string trimed = configLine;
    trim(trimed);
    size_t eqPos = trimed.find('=');
    if (eqPos == string::npos) {
        return false;
    }
    key = trimed.substr(0, eqPos);
    value = trimed.substr(eqPos + 1);
    trim(key);
    trim(value);
        return true;
}

string getStringConfigValue(string section, string key, string defultValue) {
    auto sectionIt = configData.find(section);
    if (sectionIt != configData.end()) {
        auto keyIt = sectionIt->second.find(key);
        if (keyIt != sectionIt->second.end()) {
            return keyIt->second;
        }
        else {
			return defultValue;
        }
    }
    else {
        return defultValue;
    }
}

int getIntConfigValue(string section, string key, int defultValue) {
    string valueStr = getStringConfigValue(section, key, "");
    if(valueStr.empty()){
        return defultValue;
	}
    try {
        int value = stoi(valueStr);
		return value;
    }
    catch(exception){
		return defultValue;
    }

}
double getDoubleConfigValue(string section, string key, double defultValue) {
    string valueStr = getStringConfigValue(section, key, "");
    if (valueStr.empty()) {
        return defultValue;
    }
    try {
        double value = stod(valueStr);
        return value;
    }
    catch (exception) {
        return defultValue;
    }

}
bool getBoolConfigValue(string section, string key, bool defultValue) {
    string valueStr = getStringConfigValue(section, key, "");
    if (valueStr.empty()) {
        return defultValue;
    }
    if (valueStr == "TRUE" || valueStr == "true") {
        return true;
    }
    else if (valueStr == "FALSE" || valueStr == "false") {
        return false;
    }
    else {
		return defultValue;
    }
}