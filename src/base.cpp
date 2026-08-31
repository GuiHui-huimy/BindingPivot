#include "../include/header.h"

string formatDoubleToString(double value, int precision) {
	std::ostringstream oss;
	oss << std::fixed << std::setprecision(precision) << value;
	return oss.str();
}