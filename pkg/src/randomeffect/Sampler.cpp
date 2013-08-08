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

#include "include/Sampler.h"

const double Sampler::EPSILON = 0.00000001;

Sampler::Sampler() : gwafile(NULL), reader(NULL), common_markers(auxiliary::bool_strcmp_ignore_case) {

}

Sampler::~Sampler() {
	gwafile = NULL;
	reader = NULL;

	for (common_markers_it = common_markers.begin(); common_markers_it != common_markers.end(); ++common_markers_it) {
		delete common_markers_it->second;
		free(common_markers_it->first);
	}
}

void Sampler::open_gwafile(GwaFile* gwafile) throw (SamplerException) {
	if (gwafile == NULL) {
		throw SamplerException("Sampler", "open_gwafile( GwaFile* )", __LINE__, 0, "gwafile");
	}

	try {
		close_gwafile();

		this->gwafile = gwafile;

		reader = ReaderFactory::create(gwafile->get_descriptor()->get_full_path());
		reader->set_file_name(gwafile->get_descriptor()->get_full_path());
		reader->open();
	} catch (DescriptorException& e) {
		SamplerException new_e(e);
		new_e.add_message("Sampler", "open_gwafile( GwaFile* )", __LINE__, 3, gwafile->get_descriptor()->get_full_path());
		throw new_e;
	} catch (ReaderException& e) {
		SamplerException new_e(e);
		new_e.add_message("Sampler", "open_gwafile( GwaFile* )", __LINE__, 3, gwafile->get_descriptor()->get_full_path());
		throw new_e;
	} catch (SamplerException& e) {
		e.add_message("Sampler", "open_gwafile( GwaFile* )", __LINE__, 3, gwafile->get_descriptor()->get_full_path());
		throw;
	}
}

void Sampler::close_gwafile() throw (SamplerException) {
	try {
		if (reader != NULL) {
			reader->close();
			delete reader;
			reader = NULL;
		}
	} catch (ReaderException &e) {
		SamplerException new_e(e);
		new_e.add_message("Sampler", "close_gwafile()", __LINE__, 4, (gwafile != NULL) ? gwafile->get_descriptor()->get_full_path() : "NULL");
		throw new_e;
	}

	gwafile = NULL;
}

void Sampler::process_gwafile_header() throw (SamplerException) {
	Descriptor* descriptor = NULL;

	char header_separator = '\0';
	char* header = NULL;
	char* token = NULL;
	int column_position = 0;
	const char* column_name = NULL;
	bool regions_append = false;

	if (gwafile == NULL) {
		return;
	}

	try {
		descriptor = gwafile->get_descriptor();
		header_separator = gwafile->get_header_separator();

		if (reader->read_line() <= 0) {
			throw SamplerException("Sampler", "process_gwafile_header()", __LINE__, 5, 1, gwafile->get_descriptor()->get_name());
		}

		header = *(reader->line);

		total_columns = numeric_limits<int>::min();
		marker_column_pos = numeric_limits<int>::min();
		maf_column_pos = numeric_limits<int>::min();
		oevar_imp_column_pos = numeric_limits<int>::min();

		token = auxiliary::strtok(&header, header_separator);
		while (token != NULL) {
			column_name = descriptor->get_default_column(token, gwafile->is_case_sensitive());
			if (column_name != NULL) {
				if (strcmp(column_name, Descriptor::MARKER) == 0) {
					marker_column_pos = column_position;
				} else if (strcmp(column_name, Descriptor::FREQLABEL) == 0) {
					maf_column_pos = column_position;
				} else if (strcmp(column_name, Descriptor::OEVAR_IMP) == 0) {
					oevar_imp_column_pos = column_position;
				}
			}
			token = auxiliary::strtok(&header, header_separator);
			++column_position;
		}

		total_columns = column_position;

		if (marker_column_pos < 0) {
			throw SamplerException("Sampler", "process_gwafile_header()", __LINE__, 7, ((column_name = descriptor->get_column(Descriptor::MARKER)) != NULL) ? column_name : Descriptor::MARKER, gwafile->get_descriptor()->get_name());
		}

		if (maf_column_pos < 0) {
			throw SamplerException("Sampler", "process_gwafile_header()", __LINE__, 7, ((column_name = descriptor->get_column(Descriptor::FREQLABEL)) != NULL) ? column_name : Descriptor::FREQLABEL, gwafile->get_descriptor()->get_name());
		}

		if (oevar_imp_column_pos < 0) {
			throw SamplerException("Sampler", "process_gwafile_header()", __LINE__, 7, ((column_name = descriptor->get_column(Descriptor::OEVAR_IMP)) != NULL) ? column_name : Descriptor::OEVAR_IMP, gwafile->get_descriptor()->get_name());
		}
	} catch (ReaderException &e) {
		SamplerException new_e(e);
		new_e.add_message("Sampler", "process_header_without_map()", __LINE__, 6, gwafile->get_descriptor()->get_name());
		throw new_e;
	} catch (DescriptorException &e) {
		SamplerException new_e(e);
		new_e.add_message("Sampler", "process_header_without_map()", __LINE__, 6, gwafile->get_descriptor()->get_name());
		throw new_e;
	}
}

void Sampler::initialize_common_markers() throw (SamplerException) {
	Descriptor* descriptor = NULL;
	const char* file_name = NULL;
	char data_separator = '\0';
	double oevar_imp_threshold = 0.0;
	const char* missing_value = NULL;
	const char* column_name = NULL;

	char* line = NULL;
	int line_length = 0;
	unsigned int line_number = 2u;

	char* token = NULL;
	char* end_ptr = NULL;

	int column_position = 0;

	char* marker_token = NULL;
	char* maf_token = NULL;
	char* oevar_imp_token = NULL;

	double maf_value = 0.0;
	double oevr_imp_value = 0.0;

	char* marker = NULL;
	vector<entry>* maf_values = NULL;

	entry entry_temp;

	if (gwafile == NULL) {
		return;
	}

	for (common_markers_it = common_markers.begin(); common_markers_it != common_markers.end(); ++common_markers_it) {
		delete common_markers_it->second;
		free(common_markers_it->first);
	}

	try {
		descriptor = gwafile->get_descriptor();
		file_name = descriptor->get_name();
		data_separator = gwafile->get_data_separator();
		missing_value = descriptor->get_property(Descriptor::MISSING);
		oevar_imp_threshold = descriptor->get_threshold(Descriptor::IMP)->front();

		while ((line_length = reader->read_line()) > 0) {
			line = *(reader->line);

			column_position = 0;
			marker_token = NULL;
			maf_token = NULL;
			oevar_imp_token = NULL;
			token = auxiliary::strtok(&line, data_separator);
			while (token != NULL) {
				if (column_position == marker_column_pos) {
					marker_token = token;
				} else if (column_position == maf_column_pos) {
					maf_token = token;
				} else if (column_position == oevar_imp_column_pos) {
					oevar_imp_token = token;
				}
				token = auxiliary::strtok(&line, data_separator);
				++column_position;
			}

			if (column_position < total_columns) {
				throw SamplerException("Sampler", "initialize_common_markers()", __LINE__, 8, line_number, descriptor->get_name(), column_position, total_columns);
			} else if (column_position > total_columns) {
				throw SamplerException("Sampler", "initialize_common_markers()", __LINE__, 9, line_number, descriptor->get_name(), column_position, total_columns);
			}

			if (strcmp(missing_value, oevar_imp_token) == 0) {
				++line_number;
				continue;
			}

			oevr_imp_value = strtod(oevar_imp_token, &end_ptr);
			if (*end_ptr != '\0') {
				throw SamplerException("Sampler", "initialize_common_markers()",  __LINE__, 10, oevar_imp_token, ((column_name = descriptor->get_column(Descriptor::OEVAR_IMP)) != NULL) ? column_name : Descriptor::OEVAR_IMP, line_number);
			}

			if (isnan(oevr_imp_value)) {
				++line_number;
				continue;
			}

			if (auxiliary::fcmp(oevr_imp_value, oevar_imp_threshold, EPSILON) < 0) {
				++line_number;
				continue;
			}

			if (strcmp(missing_value, maf_token) == 0) {
				++line_number;
				continue;
			}

			maf_value = strtod(maf_token, &end_ptr);
			if (*end_ptr != '\0') {
				throw SamplerException("Sampler", "initialize_common_markers()",  __LINE__, 10, maf_token, ((column_name = descriptor->get_column(Descriptor::FREQLABEL)) != NULL) ? column_name : Descriptor::FREQLABEL, line_number);
			}

			if (isnan(maf_value)) {
				++line_number;
				continue;
			}

			common_markers_it = common_markers.find(marker_token);
			if (common_markers_it != common_markers.end()) {
				throw SamplerException("Sampler", "initialize_common_markers()", __LINE__, 15, marker_token, descriptor->get_name());
			} else {
				marker = (char*)malloc((strlen(marker_token) + 1) * sizeof(char));
				if (marker == NULL) {
					throw SamplerException("Sampler", "initialize_common_markers()", __LINE__, 2, (strlen(marker_token) + 1) * sizeof(char));
				}
				strcpy(marker, marker_token);

				maf_values = new vector<entry>();

				common_markers.insert(pair<char*, vector<entry>* >(marker, maf_values));
			}

			if (auxiliary::fcmp(maf_value, 0.5, EPSILON) > 0) {
				entry_temp.maf = 1.0 - maf_value;
			} else {
				entry_temp.maf = maf_value;
			}
			entry_temp.line_id = line_number;

			maf_values->push_back(entry_temp);

			++line_number;
		}

		if (line_length == 0) {
			throw SamplerException("Sampler", "initialize_common_markers()", __LINE__, 13, line_number, descriptor->get_name());
		}

	} catch (DescriptorException &e) {
		SamplerException new_e(e);
		new_e.add_message("Sampler", "initialize_common_markers()", __LINE__, 14, gwafile->get_descriptor()->get_name());
		throw new_e;
	} catch (ReaderException &e) {
		SamplerException new_e(e);
		new_e.add_message("Sampler", "initialize_common_markers()", __LINE__, 14, gwafile->get_descriptor()->get_name());
		throw new_e;
	} catch (SamplerException &e) {
		e.add_message("Sampler", "initialize_common_markers()", __LINE__, 14, gwafile->get_descriptor()->get_name());
		throw;
	}
}

void Sampler::add_common_markers() throw (SamplerException) {
	Descriptor* descriptor = NULL;
	const char* file_name = NULL;
	char data_separator = '\0';
	double oevar_imp_threshold = 0.0;
	const char* missing_value = NULL;
	const char* column_name = NULL;

	char* line = NULL;
	int line_length = 0;
	unsigned int line_number = 2u;

	char* token = NULL;
	char* end_ptr = NULL;

	int column_position = 0;

	char* marker_token = NULL;
	char* maf_token = NULL;
	char* oevar_imp_token = NULL;

	double maf_value = 0.0;
	double oevr_imp_value = 0.0;

	entry entry_temp;

	if (gwafile == NULL) {
		return;
	}

	try {
		descriptor = gwafile->get_descriptor();
		file_name = descriptor->get_name();
		data_separator = gwafile->get_data_separator();
		missing_value = descriptor->get_property(Descriptor::MISSING);
		oevar_imp_threshold = descriptor->get_threshold(Descriptor::IMP)->front();

		while ((line_length = reader->read_line()) > 0) {
			line = *(reader->line);

			column_position = 0;
			marker_token = NULL;
			maf_token = NULL;
			oevar_imp_token = NULL;
			token = auxiliary::strtok(&line, data_separator);
			while (token != NULL) {
				if (column_position == marker_column_pos) {
					marker_token = token;
				} else if (column_position == maf_column_pos) {
					maf_token = token;
				} else if (column_position == oevar_imp_column_pos) {
					oevar_imp_token = token;
				}
				token = auxiliary::strtok(&line, data_separator);
				++column_position;
			}

			if (column_position < total_columns) {
				throw SamplerException("Sampler", "add_common_markers()", __LINE__, 8, line_number, descriptor->get_name(), column_position, total_columns);
			} else if (column_position > total_columns) {
				throw SamplerException("Sampler", "add_common_markers()", __LINE__, 9, line_number, descriptor->get_name(), column_position, total_columns);
			}

			if (strcmp(missing_value, oevar_imp_token) == 0) {
				++line_number;
				continue;
			}

			oevr_imp_value = strtod(oevar_imp_token, &end_ptr);
			if (*end_ptr != '\0') {
				throw SamplerException("Sampler", "add_common_markers()",  __LINE__, 10, oevar_imp_token, ((column_name = descriptor->get_column(Descriptor::OEVAR_IMP)) != NULL) ? column_name : Descriptor::OEVAR_IMP, line_number);
			}

			if (isnan(oevr_imp_value)) {
				++line_number;
				continue;
			}

			if (auxiliary::fcmp(oevr_imp_value, oevar_imp_threshold, EPSILON) < 0) {
				++line_number;
				continue;
			}

			if (strcmp(missing_value, maf_token) == 0) {
				++line_number;
				continue;
			}

			maf_value = strtod(maf_token, &end_ptr);
			if (*end_ptr != '\0') {
				throw SamplerException("Sampler", "add_common_markers()",  __LINE__, 10, maf_token, ((column_name = descriptor->get_column(Descriptor::FREQLABEL)) != NULL) ? column_name : Descriptor::FREQLABEL, line_number);
			}

			if (isnan(maf_value)) {
				++line_number;
				continue;
			}

			common_markers_it = common_markers.find(marker_token);
			if (common_markers_it != common_markers.end()) {
				if (auxiliary::fcmp(maf_value, 0.5, EPSILON) > 0) {
					entry_temp.maf = 1.0 - maf_value;
				} else {
					entry_temp.maf = maf_value;
				}
				entry_temp.line_id = line_number;
				common_markers_it->second->push_back(entry_temp);
			}

			++line_number;
		}

		if (line_length == 0) {
			throw SamplerException("Sampler", "add_common_markers()", __LINE__, 13, line_number, descriptor->get_name());
		}

	} catch (DescriptorException &e) {
		SamplerException new_e(e);
		new_e.add_message("Sampler", "add_common_markers()", __LINE__, 14, gwafile->get_descriptor()->get_name());
		throw new_e;
	} catch (ReaderException &e) {
		SamplerException new_e(e);
		new_e.add_message("Sampler", "add_common_markers()", __LINE__, 14, gwafile->get_descriptor()->get_name());
		throw new_e;
	} catch (SamplerException &e) {
		e.add_message("Sampler", "add_common_markers()", __LINE__, 14, gwafile->get_descriptor()->get_name());
		throw;
	}
}

void Sampler::get_common_entries(vector<GwaFile*>& gwa_files) throw (SamplerException) {
	GwaFile* gwa_file = NULL;
	vector<GwaFile*>::iterator gwa_files_it;

	vector<entry>* entries = NULL;

	gwa_files_it = gwa_files.begin();
	if (gwa_files_it != gwa_files.end()) {
		open_gwafile(*gwa_files_it);
		process_gwafile_header();
		initialize_common_markers();
		close_gwafile();
		++gwa_files_it;

		while (gwa_files_it != gwa_files.end()) {
			open_gwafile(*gwa_files_it);
			process_gwafile_header();
			add_common_markers();
			close_gwafile();
			++gwa_files_it;
		}

		common_markers_it = common_markers.begin();
		while (common_markers_it != common_markers.end()) {
			entries = common_markers_it->second;
			if (entries->size() != gwa_files.size()) {
				free(common_markers_it->first);
				delete entries;
				common_markers.erase(common_markers_it++);
			} else {
				++common_markers_it;
			}
		}
	}
}

void Sampler::filter_common_markers(map<char*, vector<entry>*, bool(*)(const char*, const char*) >& filtered_common_markers, double max_maf, bool strict) {
	char* marker = NULL;
	vector<entry>* entries = NULL;
	vector<entry>::iterator entries_it;

	double min_maf = numeric_limits<double>::max();

	if (strict) {
		common_markers_it = common_markers.begin();
		while (common_markers_it != common_markers.end()) {
			marker = common_markers_it->first;
			entries = common_markers_it->second;
			min_maf = numeric_limits<double>::max();

			for (entries_it = entries->begin(); entries_it != entries->end(); ++entries_it) {
				if (auxiliary::fcmp(entries_it->maf, min_maf, EPSILON) < 0) {
					min_maf = entries_it->maf;
				}
			}

			if (auxiliary::fcmp(min_maf, max_maf, EPSILON) < 0) {
				filtered_common_markers.insert(pair<char*, vector<entry>* >(marker, entries));
				common_markers.erase(common_markers_it++);
			} else {
				++common_markers_it;
			}
		}
	} else {
		common_markers_it = common_markers.begin();
		while (common_markers_it != common_markers.end()) {
			marker = common_markers_it->first;
			entries = common_markers_it->second;
			min_maf = numeric_limits<double>::max();

			for (entries_it = entries->begin(); entries_it != entries->end(); ++entries_it) {
				if (auxiliary::fcmp(entries_it->maf, min_maf, EPSILON) < 0) {
					min_maf = entries_it->maf;
				}
			}

			if (auxiliary::fcmp(min_maf, max_maf, EPSILON) <= 0) {
				filtered_common_markers.insert(pair<char*, vector<entry>* >(marker, entries));
				common_markers.erase(common_markers_it++);
			} else {
				++common_markers_it;
			}
		}
	}
}

void Sampler::get_random_markers(map<char*, vector<entry>*, bool(*)(const char*, const char*) >& markers, map<char*, vector<entry>*, bool(*)(const char*, const char*) >& random_markers, int n_random) {
	int i = 0;
	int n_markers = 0;
	map<char*, vector<entry>*, bool(*)(const char*, const char*) >::iterator markers_it;
	vector<bool> random;

	n_markers = markers.size();

	if (n_random > n_markers) {
		markers_it = markers.begin();
		while (markers_it != markers.end()) {
			random_markers.insert(pair<char*, vector<entry>* >(markers_it->first, markers_it->second));
			++markers_it;
		}
	} else {
		while (i < n_random) {
			random.push_back(true);
			++i;
		}

		while (i < n_markers) {
			random.push_back(false);
			++i;
		}

		random_shuffle(random.begin(), random.end());

		i = 0;
		markers_it = markers.begin();
		while (markers_it != markers.end()) {
			if (random.at(i)) {
				random_markers.insert(pair<char*, vector<entry>* >(markers_it->first, markers_it->second));
			}
			++markers_it;
			++i;
		}
	}
}

void Sampler::get_line_ids(map<char*, vector<entry>*, bool(*)(const char*, const char*) >& markers, vector<unsigned int>& line_ids, int gwafile_id) {
	map<char*, vector<entry>*, bool(*)(const char*, const char*) >::iterator markers_it;
	vector<entry>* entries = NULL;

	for (markers_it = markers.begin(); markers_it != markers.end(); ++markers_it) {
		line_ids.push_back(markers_it->second->at(gwafile_id).line_id);
	}

	sort(line_ids.begin(), line_ids.end(), std::greater<unsigned int>());
}

void Sampler::write_lines(GwaFile* gwafile, vector<unsigned int>& line_ids, const char* output_prefix) throw (SamplerException) {
	Reader* reader = NULL;

	ofstream ofile_stream;
	char* output_file_name = NULL;

	char* line = NULL;
	int line_length = 0;
	unsigned int line_number = 1u;

	if (line_ids.size() <= 0) {
		return;
	}

	try {
		auxiliary::transform_file_name(&output_file_name, output_prefix, gwafile->get_descriptor()->get_name(), NULL, true);
		if (output_file_name == NULL) {
			throw SamplerException("Sampler", "write_lines()", __LINE__, 16);
		}

		ofile_stream.exceptions(ios_base::failbit | ios_base::badbit);

		try {
			ofile_stream.open(output_file_name);
		} catch (ofstream::failure &e) {
			throw SamplerException("Sampler", "write_lines()", __LINE__, 17, output_file_name);
		}

		reader = ReaderFactory::create(gwafile->get_descriptor()->get_full_path());
		reader->set_file_name(gwafile->get_descriptor()->get_full_path());
		reader->open();

		try {
			if ((line_length = reader->read_line()) <= 0) {
				throw SamplerException("Sampler", "write_lines()", __LINE__, 5, line_number, gwafile->get_descriptor()->get_name());
			}

			ofile_stream << *(reader->line) << endl;
			++line_number;

			while ((line_length = reader->read_line()) > 0) {
				while (line_number > line_ids.back()) {
					line_ids.pop_back();
				}
				if (line_number == line_ids.back()) {
					ofile_stream << *(reader->line) << endl;
					line_ids.pop_back();
				}

				++line_number;
			}

		} catch (ofstream::failure &e) {
			throw SamplerException("Sampler", "write_lines()", __LINE__, 19, output_file_name);
		}

		reader->close();

		try {
			ofile_stream.close();
		} catch (ofstream::failure &e) {
			throw SamplerException("Sampler", "write_lines()", __LINE__, 18, output_file_name);
		}

		if (line_length == 0) {
			throw SamplerException("Sampler", "write_lines()", __LINE__, 13, line_number, gwafile->get_descriptor()->get_name());
		}
	} catch (DescriptorException &e) {
		SamplerException new_e(e);
		new_e.add_message("Sampler", "write_lines()", __LINE__, 14, gwafile->get_descriptor()->get_name());
		throw new_e;
	} catch (ReaderException &e) {
		SamplerException new_e(e);
		new_e.add_message("Sampler", "write_lines()", __LINE__, 14, gwafile->get_descriptor()->get_name());
		throw new_e;
	} catch (SamplerException &e) {
		e.add_message("Sampler", "write_lines()", __LINE__, 14, gwafile->get_descriptor()->get_name());
		throw;
	}
}

void Sampler::sample_common_markers(vector<GwaFile*>& gwa_files) throw (SamplerException) {
	map<char*, vector<entry>*, bool(*)(const char*, const char*) > filtered_common_markers(auxiliary::bool_strcmp_ignore_case);
	map<char*, vector<entry>*, bool(*)(const char*, const char*) > random_markers(auxiliary::bool_strcmp_ignore_case);
	map<char*, vector<entry>*, bool(*)(const char*, const char*) >::iterator markers_it;
	vector<unsigned int> random_line_ids;

	vector<double>* maf_thresholds = NULL;
	const char* output_prefix = NULL;

	stringstream file_prefix;

	try {
		maf_thresholds = gwa_files.back()->get_descriptor()->get_threshold(Descriptor::MAF);
		output_prefix = gwa_files.back()->get_descriptor()->get_property(Descriptor::PREFIX);
	} catch (DescriptorException &e) {
		SamplerException new_e(e);
		new_e.add_message("Sampler", "sample_common_markers()", __LINE__, 23);
		throw new_e;
	}

//	All markers
	file_prefix << output_prefix << "_ALL_";
	get_random_markers(common_markers, random_markers, 1000);

	for (int j = 0; j < gwa_files.size(); ++j) {
		get_line_ids(random_markers, random_line_ids, j);
		write_lines(gwa_files.at(j), random_line_ids, file_prefix.str().c_str());
		random_line_ids.clear();
	}

	random_markers.clear();
	file_prefix.str("");

//	Markers with MAF threshold
	for (int i = 0; i < maf_thresholds->size(); ++i) {
		file_prefix << output_prefix;
		if (i > 0) {
			file_prefix  << "_FROM_" << maf_thresholds->at(i - 1);
		}
		file_prefix  << "_UPTO_" << maf_thresholds->at(i) << "_";

		filter_common_markers(filtered_common_markers, maf_thresholds->at(i), i > 0 ? false : true );
		get_random_markers(filtered_common_markers, random_markers, 1000);

		for (int j = 0; j < gwa_files.size(); ++j) {
			get_line_ids(random_markers, random_line_ids, j);
			write_lines(gwa_files.at(j), random_line_ids, file_prefix.str().c_str());
			random_line_ids.clear();
		}

		random_markers.clear();
		file_prefix.str("");

		markers_it = filtered_common_markers.begin();
		while (markers_it != filtered_common_markers.end()) {
			free(markers_it->first);
			delete markers_it->second;
			++markers_it;
		}
		filtered_common_markers.clear();
	}

//	All the rest markers
	if (common_markers.size() > 0) {
		file_prefix << output_prefix << "_FROM_" <<  maf_thresholds->back() << "_";
		get_random_markers(common_markers, random_markers, 1000);

		for (int j = 0; j < gwa_files.size(); ++j) {
			get_line_ids(random_markers, random_line_ids, j);
			write_lines(gwa_files.at(j), random_line_ids, file_prefix.str().c_str());
			random_line_ids.clear();
		}

		random_markers.clear();
		file_prefix.str("");
	}
}
