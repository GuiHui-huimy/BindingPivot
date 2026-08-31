#include "../include/header.h"
int main() {
	ErrorHeader ErrorInfo;
	string workPath = "";
	string pdbFullPath = "";
	
	ErrorInfo =loadConfig();
	if (ErrorInfo.isError) {
		cerr << "Error: " << ErrorInfo.errorMessage << endl;
	}
	ErrorInfo =IndexPage(workPath);
	if (ErrorInfo.isError) {
		cerr << "Error: " << ErrorInfo.errorMessage << endl;
	}
	while (1){
		ErrorInfo = choseList(workPath);
		if (ErrorInfo.isError) {
			cerr << "Error: " << ErrorInfo.errorMessage << endl;
			continue;
		}
	}
	return 0;
}




