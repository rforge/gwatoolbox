/*
 * Copyright © 2012 Daniel Taliun, Christian Fuchsberger and Cristian Pattaro. All rights reserved.
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

#ifndef SELECTOR_H_
#define SELECTOR_H_

#include "SelectorException.h"
#include "../../gwafile/include/GwaFile.h"

class Selector {
private:
	GwaFile* gwafile;

	Reader* reader;

	char* header_backup;

	int total_columns;
	int marker_column_pos;
	int chr_column_pos;
	int pvalue_column_pos;

public:
	Selector();
	virtual ~Selector();

	void open_gwafile(GwaFile* gwafile) throw (SelectorException);
	void close_gwafile() throw (SelectorException);

	void process_header() throw (SelectorException);
};

#endif
