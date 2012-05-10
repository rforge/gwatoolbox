/*
 * Copyright © 2011 Daniel Taliun, Christian Fuchsberger and Cristian Pattaro. All rights reserved.
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

#ifndef ANNOTATOR_H_
#define ANNOTATOR_H_

#include <map>
#include <vector>

#include "AnnotatorException.h"
#include "../../auxiliary/include/auxiliary.h"
#include "../../reader/include/TextReader.h"
#include "../../gwafile/include/GwaFile.h"
#include "../../intervaltree/include/IntervalTree.h"

class Annotator {
private:
	GwaFile* gwafile;
	const char* regions_file;

	TextReader reader;
	TextReader regions_reader;

	char* header_backup;

	int total_columns;
	int marker_column_pos;
	int chr_column_pos;
	int position_column_pos;

	int regions_file_total_columns;
	int region_name_column_pos;
	int region_chr_column_pos;
	int region_start_column_pos;
	int region_end_column_pos;

	map<char*, IntervalTree<char*>*, bool(*)(const char*, const char*)> regions_indices;
	map<char*, IntervalTree<char*>*, bool(*)(const char*, const char*)>::iterator regions_indices_it;

	void open_regions_file() throw (AnnotatorException);
	void close_regions_file() throw (AnnotatorException);
	void process_regions_file_header() throw (AnnotatorException);
	void process_regions_file_data() throw (AnnotatorException);

	void write_char_vector(ofstream &ofile_stream, vector<char*>* values, char separator) throw (ofstream::failure);

public:
	Annotator();
	virtual ~Annotator();

	void open_gwafile(GwaFile* gwafile) throw (AnnotatorException);
	void close_gwafile() throw (AnnotatorException);

	void process_header() throw (AnnotatorException);
	void index_regions() throw (AnnotatorException);
	void annotate() throw (AnnotatorException);
};

#endif
