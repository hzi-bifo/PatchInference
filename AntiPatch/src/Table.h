/*
 * Helpers.h
 *
 *  Created on: Jun 10, 2012
 *      Author: ctusche
 */

#ifndef TABLE_H
#define TABLE_H

#include <string>
#include <vector>
#include <iostream>



class Table {

private:
	std::vector< std::string > entries;
	void create( std::string sep );

public:

	std::vector< std::vector<std::string> > T;

	Table();

	Table(std::string path, std::string sep=" \t");
	void print(std::ostream &target=std::cout);
	unsigned long rows();
	std::vector< std::string > col(unsigned long long  k);
	std::vector< std::string > row(unsigned long long  k);


};


#endif /* TABLE_H */
