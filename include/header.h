// The code is part of a project that deals with atomic structures and calculations related to atomic sites.
// header.h
#ifndef HEADER_H
#define HEADER_H
#define ERROR 0
#define SUCCESS 1
//Include files
#include <iostream> 
#include <vector>
#include <regex>
#include <map>
#include <sstream>
#include <string>
#include <set>
#include <iomanip> 
#include <fstream>
using namespace std;

// Structure to hold atom information
struct AtomInfo {
	string recordName;
	string moleculeName;
	int moleculeSerial;
	string atomName; // 原子名称（如CA、N、O）
	string type;     // 原子类型（如C、O、N）
	double x, y, z; // 坐标
	AtomInfo() : recordName(""), moleculeName(""), moleculeSerial(0),
		atomName(""), type(""), x(0.0), y(0.0), z(0.0) {
	}
};

typedef struct ErrorHeader{
	bool isError = false;
	string errorMessage = "";
	void ErrorLogOut() {
		cout << "Error: " << errorMessage << endl;
	}
}ErrorHeader;


// ---Function declarations

//rWeightedSiteCalc functions
bool initArrayDect(int atomNum, std::vector<string> atomType, std::vector<double> x, std::vector<double> y, std::vector<double> z);
std::vector<double> rWeightedSiteCalc(int atomNum, std::vector<string> atomType, std::vector<double> x, std::vector<double> y, std::vector<double> z);

//pdb analysis functions
double getAtomicRadius(string atomType);
bool extractSingleLine(const string& line, struct AtomInfo& info);
string extractElementFromAtomName(const string& atomName);
bool specialLine(string line, bool isOutput);
string smartElementExtraction(const string& atomName);
/*vector<string> getLinesfromMutiString(const string& pdbData);*/
int initGetInfoFromLines(rsize_t& lineNum, string fileDir, vector<string>& Lines, vector<double>& x, vector<double>& y, vector<double>& z, vector<string>& atomType);
string createPdbLine(double x, double y, double z);


//core analysis funtions
ErrorHeader analysisLigandCoreSite(string fileDir);
ErrorHeader analysisResCoreSite(string fileDir);
ErrorHeader analysisPocketCoreSite(string fileDir);
bool isNeedMol(string MolName, int MolSer, AtomInfo infos);
ErrorHeader recogLigand(int lineNum, vector<string> Lines, vector <string>& recLigands, vector<int>& recLigandseries);

//index page functions
ErrorHeader IndexPage(string &workPath);
ErrorHeader setWorkPath(string& workPath);
ErrorHeader chosePdbAnalysis(string& workPath, string& fullPdbPath);
ErrorHeader choseList(string& workPath);


//file functions
int getConfigData(string fileDir, vector<string>& configLines);
int getPdbData(string fileDir, vector<string>& pdbLines);
int getProteinList(string fileDir, vector<string>& proteinList);
int writeCsvData(string fileDir, vector<string> initRows, vector<string> rows);
ErrorHeader writeDefaultBatchConfig(string workPath, bool& creatDefault);

//config functions
ErrorHeader loadConfig();
void trim(std::string& str);
bool isSectionHeader(string configLine, string& sectionName);
bool isCommendLine(string configLine);
bool pharseKeyValue(string configLine, string& key, string& value);
ErrorHeader readConfigData(string ConfigPath);
string getStringConfigValue(string section, string key, string defultValue);
int getIntConfigValue(string section, string key, int defultValue);
double getDoubleConfigValue(string section, string key, double defultValue);
bool getBoolConfigValue(string section, string key, bool defultValue);
string getExePath();

//batch functions
ErrorHeader ligandBatchCore(string fileDir, string ligandProteinName, string modelName, vector <string>& returnrecLigands, vector<double>& coreSite);
ErrorHeader batchProcess(string workPath);
void batchInstructions(string workPath);
ErrorHeader pocketBatchCore(string fileDir, string pocketProteinName, string modelName, vector<double>& returnCoreSite);
bool batchConfigInfo();

//base functions
string formatDoubleToString(double value, int precision = 3);

//Define constants
const set<string> VALID_ELEMENTS = {
	"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
	"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
	"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
	"Ga", "Ge", "As", "Se", "Br", "Kr", "I", "Pt", "Au", "Ag"
};
//Map we need special end set !!
const map<string,string> KNOW_NOMEAN_R = {
	{"MIX","NOMEAN"}
};
static const map<string, double>	RADUIS_TABLE = {
	{"H", 1.20}, {"C", 1.70}, {"N", 1.55}, {"O", 1.52},
	{"F", 1.47}, {"S", 1.80}, {"Cl", 1.75}, {"P", 1.80},
	{"Br", 1.85}, {"I", 1.98}, {"Fe", 1.40}, {"Zn", 1.39},
	{"Cu", 1.40}, {"Mg", 1.73}, {"Ca", 2.31}, {"Na", 2.27},
	{"K", 2.75},{"NOMEAN",1.00}
};

static vector<string> SPECIAL_RECORD = {
	"HEADER","TITLE","COMPND","SOURCE","KEYWDS","EXPDTA",
	"REMARK","DBREF","SEQRES","CRYST1","ORIGX1","ORIGX2",
	"ORIGX3","SCALE1","SCALE2","SCALE3","MASTER","MODEL",
	"CONECT"
};

extern map<string, map<string, string>> configData;


#endif