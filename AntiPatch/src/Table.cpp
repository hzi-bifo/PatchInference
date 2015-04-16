/*
 * Helpers.cpp
 *
 *  Created on: Jun 10, 2012
 *      Author: ctusche
 */

#include <unistd.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include "Table.h"

std::vector<std::string> tokenizeToString(const std::string& str,
		const std::string& delimiters) {

	std::vector<std::string> tokens;
	std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	std::string::size_type pos = str.find_first_of(delimiters, lastPos);
	std::string tmp;

	if (std::string::npos == pos || std::string::npos == lastPos) {
		tokens.push_back(str);
		return tokens;
	}

	while (std::string::npos != pos || std::string::npos != lastPos) {
		tmp = str.substr(lastPos, pos - lastPos);
		tokens.push_back(tmp);
		lastPos = str.find_first_not_of(delimiters, pos);
		pos = str.find_first_of(delimiters, lastPos);
	}

	return tokens;
}

Table::Table(std::string path, std::string sep) {
	std::ifstream in;
	in.open(path.c_str());
	std::vector<std::string> parts;
	std::vector<std::vector<std::string> > res;

	if (!in.is_open()) {
		std::string msg("Table (constructor): unable to open file \n\t");
		msg.append(path);
		throw msg;
	}

	std::string line;
	while (!in.eof()) {
		getline(in, line);
		if (line[0] != '#') {
			parts = tokenizeToString(line, sep);
			res.push_back(parts);
		}
	}
	this->T = res;
}

Table::Table() {}

void Table::print(std::ostream &target) {

	for (size_t i = 0; i < this->T.size(); i++) {
		for (size_t j = 0; j < T[i].size(); j++) {
			target << T[i][j] << "\t";
		}
		target << std::endl;
	}

}

unsigned long Table::rows() {
	return this->T.size();
}

std::vector<std::string> Table::col(unsigned long long k) {

	std::vector<std::string> res;
	for (unsigned long long i = 0; i < this->T.size(); i++) {
		if (this->T[i].size() > k) {
			res.push_back(this->T[i][k]);
		}
	}
	return res;
}

std::vector<std::string> Table::row(unsigned long long k) {
	return this->T[k];
}
