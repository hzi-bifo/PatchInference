/* 
 * File:   Residue.cpp
 * Author: l.Klesper
 * 
 * Created on 2. August 2012, 12:38
 */

#include "Residue.h"
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <math.h>

Residue::Residue(std::string pos, Point c, float pval, std::string alias,
		std::string epi) {

	this->pos = pos;
	this->pValue = pval;
	this->coord = c;
	this->alias = alias;
	this->epi = epi;
	this->patched = false;
	this->surface = false;
	this->agw = 0.0f;

	std::stringstream oss;
	oss << "[label=\"" << pos;
	oss << "\\n" << pval;
	if (epi != "")
		oss << "  " << epi << "\",shape=box,style=filled,color=\"0 1 1\"]";
	else
		oss << "\",shape=box,style=filled,color=\".5 .5 .5\"]";

	this->label = oss.str();
}

void Residue::setPValue(float value) {
	this->pValue = value;
}

void Residue::setAgW(float value) {
	this->agw = value;
}

Point Residue::getCoord() {
	return this->coord;
}

void Residue::setCoord(float x, float y, float z) {
	this->coord.setPoint(x, y, z);
}

float Residue::getPValue() {
	return this->pValue;
}

float Residue::getAgW() {
	return this->agw;
}

std::string Residue::getAlias() {
	return this->alias;
}

std::string Residue::getEpi() {
	return this->epi;
}

std::string Residue::getLabel() {
	return this->label;
}

std::string Residue::getPos() {
	return this->pos;
}

float Residue::getWeightedDistance(Residue* res) {
	float w = pValue - res->pValue, x = coord.getX() - res->coord.getX(), y =
			coord.getY() - res->coord.getY(), z = coord.getZ()
			- res->coord.getZ();
	return x * x + y * y + z * z + w * w;
}

float Residue::getDistance(Residue* res) {
	float x = coord.getX() - res->coord.getX(), y = coord.getY()
			- res->coord.getY(), z = coord.getZ() - res->coord.getZ();
	return x * x + y * y + z * z;
}

float Residue::getEuklideanDistance(Residue* res) {
	float x = coord.getX() - res->coord.getX(), y = coord.getY()
			- res->coord.getY(), z = coord.getZ() - res->coord.getZ();
	return sqrt(x * x + y * y + z * z);
}

bool Residue::isPatched() {
	return this->patched;
}

void Residue::setPatched(bool j) {
	this->patched = j;
}
