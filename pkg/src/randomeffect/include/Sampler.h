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

#ifndef SAMPLER_H_
#define SAMPLER_H_

#include <map>
#include <vector>
#include <algorithm>
#include <sstream>
#include <functional>

#include "SamplerException.h"
#include "../../gwafile/include/GwaFile.h"

using namespace std;

class Sampler {
private:
	GwaFile* gwafile;

	Reader* reader;

	int total_columns;
	int marker_column_pos;
	int maf_column_pos;
	int oevar_imp_column_pos;

	struct entry {
		double maf;
		unsigned int line_id;
	};

	map<char*, vector<entry>*, bool(*)(const char*, const char*)> common_markers;
	map<char*, vector<entry>*, bool(*)(const char*, const char*)>::iterator common_markers_it;

	void open_gwafile(GwaFile* gwafile) throw (SamplerException);
	void close_gwafile() throw (SamplerException);
	void process_gwafile_header() throw (SamplerException);
	void initialize_common_markers() throw (SamplerException);
	void add_common_markers() throw (SamplerException);
	void filter_common_markers(map<char*, vector<entry>*, bool(*)(const char*, const char*) >& filtered_common_markers, double max_maf, bool strict);
	void get_random_markers(map<char*, vector<entry>*, bool(*)(const char*, const char*) >& markers, map<char*, vector<entry>*, bool(*)(const char*, const char*) >& random_markers, int n_random);
	void get_line_ids(map<char*, vector<entry>*, bool(*)(const char*, const char*) >& markers, vector<unsigned int>& line_ids, int gwafile_id);
	void write_lines(GwaFile* gwafile, vector<unsigned int>& line_ids, const char* output_prefix) throw (SamplerException);

public:
	static const double EPSILON;

	Sampler();
	virtual ~Sampler();

	void get_common_entries(vector<GwaFile*>& gwa_files) throw (SamplerException);
	void sample_common_markers(vector<GwaFile*>& gwa_files) throw (SamplerException);
};

#endif
