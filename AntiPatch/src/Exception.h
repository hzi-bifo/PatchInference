/* 
 * File:   Exception.h
 * Author: l.Klesper
 *
 * Created on 2. August 2012, 18:46
 */

#ifndef EXCEPTION_H
#define	EXCEPTION_H

#include <sstream>

class Exception {
public:
    std::string error;

    Exception(std::string arg) : error(arg) {}

    template <class T>
    Exception(std::string arg1, T arg2) {
        std::stringstream ss;
        ss << arg1;
        ss << arg2;
        error = ss.str();
    }
};

#endif	/* EXCEPTION_H */


