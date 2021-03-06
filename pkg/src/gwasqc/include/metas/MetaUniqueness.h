/*
 * Copyright � 2011 Daniel Taliun, Christian Fuchsberger and Cristian Pattaro. All rights reserved.
 *
 * This file is part of GWAtoolbox.
 *
 * GWAtoolbox is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GWAtoolbox is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GWAtoolbox.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef METAUNIQUENESS_H_
#define METAUNIQUENESS_H_

#include <vector>
#include <iostream>
#include <iomanip>
#include <set>

#include "../../../auxiliary/include/auxiliary.h"
#include "Meta.h"

using namespace auxiliary;

class MetaUniqueness: public Meta {
private:
	int n;
	bool na_value;
	char** data;
	char** new_data;
	char* new_value;
	int current_heap_size;
	vector<char*> duplicates;
	vector<char*>::iterator duplicates_it;

public:
	MetaUniqueness(unsigned int heap_size = Meta::HEAP_SIZE) throw (MetaException);
	virtual ~MetaUniqueness();
	void put(char* value) throw (MetaException);
	void finalize() throw (MetaException);
	bool is_na();
	void print(ostream& stream);
	void print_html(ostream& stream, char path_separator);
	double get_memory_usage();
};

#endif
