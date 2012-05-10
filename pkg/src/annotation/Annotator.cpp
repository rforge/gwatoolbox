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

#include "include/Annotator.h"

Annotator::Annotator() : gwafile(NULL), regions_file(NULL),
	header_backup(NULL),
	total_columns(numeric_limits<int>::min()),
	marker_column_pos(numeric_limits<int>::min()),
	chr_column_pos(numeric_limits<int>::min()),
	position_column_pos(numeric_limits<int>::min()),
	regions_file_total_columns(numeric_limits<int>::min()),
	region_name_column_pos(numeric_limits<int>::min()),
	region_chr_column_pos(numeric_limits<int>::min()),
	region_start_column_pos(numeric_limits<int>::min()),
	region_end_column_pos(numeric_limits<int>::min()),
	regions_indices(auxiliary::bool_strcmp) {

}

Annotator::~Annotator() {
	regions_indices_it = regions_indices.begin();
	while (regions_indices_it != regions_indices.end()) {
		free(regions_indices_it->first);
		delete regions_indices_it->second;
		regions_indices_it++;
	}
	regions_indices.clear();

	if (header_backup != NULL) {
		free(header_backup);
		header_backup = NULL;
	}

	gwafile = NULL;
	regions_file = NULL;
}

void Annotator::open_gwafile(GwaFile* gwafile) throw (AnnotatorException) {
	if (gwafile == NULL) {
		throw AnnotatorException("Annotator", "open_gwafile( GwaFile* )", __LINE__, 0, "gwafile");
	}

	try {
		close_gwafile();
		this->gwafile = gwafile;
		reader.set_file_name(gwafile->get_descriptor()->get_full_path());
		reader.open();
	} catch (ReaderException& e) {
		AnnotatorException new_e(e);
		new_e.add_message("Annotator", "open_gwafile( GwaFile* )", __LINE__, 3, gwafile->get_descriptor()->get_full_path());
		throw new_e;
	} catch (AnnotatorException& e) {
		e.add_message("Annotator", "open_gwafile( GwaFile* )", __LINE__, 3, gwafile->get_descriptor()->get_full_path());
		throw;
	}
}

void Annotator::close_gwafile() throw (AnnotatorException) {
	try {
		reader.close();
	} catch (ReaderException &e) {
		AnnotatorException new_e(e);
		new_e.add_message("Annotator", "close_gwafile()", __LINE__, 4, (gwafile != NULL) ? gwafile->get_descriptor()->get_full_path() : "NULL");
		throw new_e;
	}

	regions_indices_it = regions_indices.begin();
	while (regions_indices_it != regions_indices.end()) {
		free(regions_indices_it->first);
		delete regions_indices_it->second;
		regions_indices_it++;
	}
	regions_indices.clear();

	if (header_backup != NULL) {
		free(header_backup);
		header_backup = NULL;
	}

	gwafile = NULL;
}

void Annotator::process_header() throw (AnnotatorException) {
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
		regions_append = gwafile->is_regions_append_on();

		if (reader.read_line() <= 0) {
			throw AnnotatorException("Annotator", "process_header()", __LINE__, 5, 1, gwafile->get_descriptor()->get_name());
		}

		header = *reader.line;

		if (regions_append) {
			header_backup = (char*)malloc((strlen(header) + 1u) * sizeof(char));
			if (header_backup == NULL) {
				throw AnnotatorException("Annotator", "process_header", __LINE__, 2, ((strlen(header) + 1u) * sizeof(char)));
			}
			strcpy(header_backup, header);
		}


		total_columns = numeric_limits<int>::min();
		marker_column_pos = numeric_limits<int>::min();
		chr_column_pos = numeric_limits<int>::min();
		position_column_pos = numeric_limits<int>::min();

		token = auxiliary::strtok(&header, header_separator);
		while (token != NULL) {
			column_name = descriptor->get_default_column(token, gwafile->is_case_sensitive());
			if (column_name != NULL) {
				if (strcmp(column_name, Descriptor::MARKER) == 0) {
					marker_column_pos = column_position;
				} else if (strcmp(column_name, Descriptor::CHR) == 0) {
					chr_column_pos = column_position;
				} else if (strcmp(column_name, Descriptor::POSITION) == 0) {
					position_column_pos = column_position;
				}
			}
			token = auxiliary::strtok(&header, header_separator);
			++column_position;
		}

		total_columns = column_position;

		if (marker_column_pos < 0) {
			throw AnnotatorException("Annotator", "process_header()", __LINE__, 7, ((column_name = descriptor->get_column(Descriptor::MARKER)) != NULL) ? column_name : Descriptor::MARKER, gwafile->get_descriptor()->get_name());
		}

		if (chr_column_pos < 0) {
			throw AnnotatorException("Annotator", "process_header()", __LINE__, 7, ((column_name = descriptor->get_column(Descriptor::CHR)) != NULL) ? column_name : Descriptor::CHR, gwafile->get_descriptor()->get_name());
		}

		if (position_column_pos < 0) {
			throw AnnotatorException("Annotator", "process_header()", __LINE__, 7, ((column_name = descriptor->get_column(Descriptor::POSITION)) != NULL) ? column_name : Descriptor::POSITION, gwafile->get_descriptor()->get_name());
		}
	} catch (ReaderException &e) {
		AnnotatorException new_e(e);
		new_e.add_message("Annotator", "process_header()", __LINE__, 6, gwafile->get_descriptor()->get_name());
		throw new_e;
	} catch (DescriptorException &e) {
		AnnotatorException new_e(e);
		new_e.add_message("Annotator", "process_header()", __LINE__, 6, gwafile->get_descriptor()->get_name());
		throw new_e;
	}
}

void Annotator::index_regions() throw (AnnotatorException) {
	if (gwafile == NULL) {
		return;
	}

	regions_indices_it = regions_indices.begin();
	while (regions_indices_it != regions_indices.end()) {
		free(regions_indices_it->first);
		delete regions_indices_it->second;
		regions_indices_it++;
	}
	regions_indices.clear();

	try {
		open_regions_file();
		process_regions_file_header();
		process_regions_file_data();
		close_regions_file();
	} catch (AnnotatorException &e) {
		e.add_message("Annotator", "index_regions()", __LINE__, 15);
		throw;
	}
}

void Annotator::annotate() throw (AnnotatorException) {
	Descriptor* descriptor = NULL;

	ofstream ofile_stream;

	const char* output_prefix = NULL;
	const char* file_name = NULL;
	char* output_file_name = NULL;

	char* line = NULL;
	int line_length = 0;
	unsigned int line_number = 2u;

	char header_separator = '\0';
	char data_separator = '\0';
	char region_separator = '\0';
	bool regions_append = false;

	char* token = NULL;
	char* end_ptr = NULL;

	int column_position = 0;

	const char* column_name = NULL;

	char* marker_token = NULL;
	char* chr_token = NULL;
	char* position_token = NULL;

	int position = 0;

	vector<double>* deviation = NULL;
	int deviation_value = 0;

	IntervalTree<char*>* genes_index = NULL;
	IntervalTree<char*> genes_index_subset;

	map<int, vector<char*>*> annotated_genes;
	map<int, vector<char*>*>::iterator annotated_genes_it;

	if (gwafile == NULL) {
		return;
	}

	try {
		descriptor = gwafile->get_descriptor();
		output_prefix = descriptor->get_property(Descriptor::PREFIX);
		file_name = descriptor->get_name();
		header_separator = gwafile->get_header_separator();
		data_separator = gwafile->get_data_separator();
		deviation = descriptor->get_threshold(Descriptor::REGIONS_DEVIATION);
		regions_append = gwafile->is_regions_append_on();

		auxiliary::transform_file_name(&output_file_name, output_prefix, file_name, NULL, true);
		if (output_file_name == NULL) {
			throw AnnotatorException("Annotator", "annotate()", __LINE__, 16);
		}

		if (data_separator == ',') {
			region_separator = ';';
		} else {
			region_separator = ',';
		}

		ofile_stream.exceptions(ios_base::failbit | ios_base::badbit);

		try {
			ofile_stream.open(output_file_name);
		} catch (ofstream::failure &e) {
			throw AnnotatorException("Annotator", "annotate()", __LINE__, 17, output_file_name);
		}

		try {
			if (regions_append) {
				ofile_stream << header_backup;
			} else {
				ofile_stream << ((column_name = descriptor->get_column(Descriptor::MARKER)) != NULL ? column_name : Descriptor::MARKER);
				ofile_stream << header_separator;
				ofile_stream << ((column_name = descriptor->get_column(Descriptor::CHR)) != NULL ? column_name : Descriptor::CHR);
				ofile_stream << header_separator;
				ofile_stream << ((column_name = descriptor->get_column(Descriptor::POSITION)) != NULL ? column_name : Descriptor::POSITION);
			}
			for (unsigned int i = 0u; i < deviation->size(); ++i) {
				deviation_value = (int)deviation->at(i);
				if (deviation_value != 0) {
					ofile_stream << header_separator << "+/-" << deviation_value;
				} else {
					ofile_stream << header_separator << "IN";
				}
			}
			ofile_stream << endl;

			while ((line_length = reader.read_line()) > 0) {
				line = *reader.line;

				if (regions_append) {
					ofile_stream << line;
				}

				column_position = 0;
				marker_token = NULL;
				chr_token = NULL;
				position_token = NULL;
				token = auxiliary::strtok(&line, data_separator);
				while (token != NULL) {
					if (column_position == marker_column_pos) {
						marker_token = token;
					} else if (column_position == chr_column_pos) {
						chr_token = token;
					} else if (column_position == position_column_pos) {
						position_token = token;
					}
					token = auxiliary::strtok(&line, data_separator);
					++column_position;
				}

				if (column_position < total_columns) {
					throw AnnotatorException("Annotator", "annotate()", __LINE__, 8, line_number, gwafile->get_descriptor()->get_name(), column_position, total_columns);
				} else if (column_position > total_columns) {
					throw AnnotatorException("Annotator", "annotate()", __LINE__, 9, line_number, gwafile->get_descriptor()->get_name(), column_position, total_columns);
				}

				position = (int)strtol(position_token, &end_ptr, 10);
				if (*end_ptr != '\0') {
					throw AnnotatorException("Annotator", "annotate()",  __LINE__, 10, position, ((column_name = descriptor->get_column(Descriptor::POSITION)) != NULL) ? column_name : Descriptor::POSITION, line_number);
				}

				if (position < 0) {
					throw AnnotatorException("Annotator", "annotate()",  __LINE__, 11, ((column_name = descriptor->get_column(Descriptor::POSITION)) != NULL) ? column_name : Descriptor::POSITION, position, line_number);
				}

				if (!regions_append) {
					ofile_stream << marker_token << data_separator << chr_token << data_separator << position_token;
				}

				regions_indices_it = regions_indices.find(chr_token);
				if (regions_indices_it != regions_indices.end()) {
					genes_index = regions_indices_it->second;

					deviation_value = (int)deviation->back();
					genes_index->get_intersecting_intervals(position - deviation_value, position + deviation_value, deviation_value, genes_index_subset);
					for (int i = deviation->size() - 2; i >= 0; --i) {
						deviation_value = (int)deviation->at(i);
						genes_index_subset.mark_intersecting_intervals(position - deviation_value, position + deviation_value, deviation_value);
					}

					genes_index_subset.get_marked_values(annotated_genes);

					if (annotated_genes.size() > 0) {
						for (unsigned int i = 0u; i < deviation->size(); ++i) {
							ofile_stream << data_separator;
							annotated_genes_it = annotated_genes.find(((int)deviation->at(i)));
							if (annotated_genes_it != annotated_genes.end()) {
								write_char_vector(ofile_stream, annotated_genes_it->second, ',');
							}
						}
						ofile_stream << endl;
					} else {
						for (unsigned int i = 0u; i < deviation->size(); ++i) {
							ofile_stream << data_separator;
						}
						ofile_stream << endl;
					}

					annotated_genes_it = annotated_genes.begin();
					while(annotated_genes_it != annotated_genes.end()) {
						delete annotated_genes_it->second;
						annotated_genes_it++;
					}
					annotated_genes.clear();

					genes_index_subset.clear();
				} else {
					for (unsigned int i = 0u; i < deviation->size(); ++i) {
						ofile_stream << data_separator;
					}
					ofile_stream << endl;
				}

				++line_number;
			}
		} catch (ofstream::failure &e) {
			throw AnnotatorException("Annotator", "annotate()", __LINE__, 19, output_file_name);
		}

		try {
			ofile_stream.close();
		} catch (ofstream::failure &e) {
			throw AnnotatorException("Annotator", "annotate()", __LINE__, 18, output_file_name);
		}

		if (line_length == 0) {
			throw AnnotatorException("Annotator", "annotate()", __LINE__, 13, line_number, gwafile->get_descriptor()->get_name());
		}
	} catch (DescriptorException &e) {
		AnnotatorException new_e(e);
		new_e.add_message("Annotator", "annotate()", __LINE__, 14, gwafile->get_descriptor()->get_name());
		throw new_e;
	} catch (ReaderException &e) {
		AnnotatorException new_e(e);
		new_e.add_message("Annotator", "annotate()", __LINE__, 14, gwafile->get_descriptor()->get_name());
		throw new_e;
	} catch (AnnotatorException &e) {
		e.add_message("Annotator", "annotate()", __LINE__, 14, gwafile->get_descriptor()->get_name());
		throw;
	}
}

void Annotator::open_regions_file() throw (AnnotatorException) {
	if (gwafile == NULL) {
		return;
	}

	try {
		close_regions_file();

		regions_file = gwafile->get_descriptor()->get_property(Descriptor::REGIONS_FILE);
		if (regions_file == NULL) {
			return;
		}

		regions_reader.set_file_name(regions_file);
		regions_reader.open();
	} catch (ReaderException& e) {
		AnnotatorException new_e(e);
		new_e.add_message("Annotator", "open_regions_file()", __LINE__, 3, (regions_file != NULL) ? regions_file : "NULL");
		throw new_e;
	}  catch (DescriptorException& e) {
		AnnotatorException new_e(e);
		new_e.add_message("Annotator", "open_regions_file()", __LINE__, 3, (regions_file != NULL) ? regions_file : "NULL");
		throw new_e;
	} catch (AnnotatorException& e) {
		e.add_message("Annotator", "open_regions_file()", __LINE__, 3, (regions_file != NULL) ? regions_file : "NULL");
		throw;
	}
}

void Annotator::close_regions_file() throw (AnnotatorException) {
	try {
		regions_reader.close();
	} catch (ReaderException &e) {
		AnnotatorException new_e(e);
		new_e.add_message("Annotator", "close_regions_file()", __LINE__, 4, (regions_file != NULL) ? regions_file : "NULL");
		throw new_e;
	}

	regions_file = NULL;
}

void Annotator::process_regions_file_header() throw (AnnotatorException) {
	Descriptor* descriptor = NULL;
	char header_separator = '\0';
	char* header = NULL;
	char* token = NULL;

	int column_position = 0;

	const char* region_name_column = NULL;
	const char* region_chr_column = NULL;
	const char* region_start_column = NULL;
	const char* region_end_column = NULL;

	if  ((gwafile == NULL) || (regions_file == NULL)) {
		return;
	}

	try {
		descriptor = gwafile->get_descriptor();
		header_separator = gwafile->get_regions_file_header_separator();

		if (regions_reader.read_line() <= 0) {
			throw AnnotatorException("Annotator", "process_regions_file_header()", __LINE__, 5, 1, regions_file);
		}

		header = *regions_reader.line;

		region_name_column = descriptor->get_property(Descriptor::REGION_NAME);
		region_chr_column = descriptor->get_property(Descriptor::REGION_CHR);
		region_start_column = descriptor->get_property(Descriptor::REGION_START);
		region_end_column = descriptor->get_property(Descriptor::REGION_END);

		regions_file_total_columns = numeric_limits<int>::min();
		region_name_column_pos = numeric_limits<int>::min();
		region_chr_column_pos = numeric_limits<int>::min();
		region_start_column_pos = numeric_limits<int>::min();
		region_end_column_pos = numeric_limits<int>::min();

		token = auxiliary::strtok(&header, header_separator);
		if (gwafile->is_case_sensitive()) {
			while (token != NULL) {
				if (strcmp(token, region_name_column) == 0) {
					region_name_column_pos = column_position;
				} else if (strcmp(token, region_chr_column) == 0) {
					region_chr_column_pos = column_position;
				} else if (strcmp(token, region_start_column) == 0) {
					region_start_column_pos = column_position;
				} else if (strcmp(token, region_end_column) == 0) {
					region_end_column_pos = column_position;
				}
				token = auxiliary::strtok(&header, header_separator);
				++column_position;
			}
		} else {
			while (token != NULL) {
				if (auxiliary::strcmp_ignore_case(token, region_name_column) == 0) {
					region_name_column_pos = column_position;
				} else if (auxiliary::strcmp_ignore_case(token, region_chr_column) == 0) {
					region_chr_column_pos = column_position;
				} else if (auxiliary::strcmp_ignore_case(token, region_start_column) == 0) {
					region_start_column_pos = column_position;
				} else if (auxiliary::strcmp_ignore_case(token, region_end_column) == 0) {
					region_end_column_pos = column_position;
				}
				token = auxiliary::strtok(&header, header_separator);
				++column_position;
			}
		}

		regions_file_total_columns = column_position;

		if (region_name_column_pos < 0) {
			throw AnnotatorException("Annotator", "process_regions_file_header()", __LINE__, 7, (region_name_column != NULL) ? region_name_column : "NULL", regions_file);
		}

		if (region_chr_column_pos < 0) {
			throw AnnotatorException("Annotator", "process_regions_file_header()", __LINE__, 7, (region_chr_column != NULL) ? region_chr_column : "NULL", regions_file);
		}

		if (region_start_column_pos < 0) {
			throw AnnotatorException("Annotator", "process_regions_file_header()", __LINE__, 7, (region_start_column != NULL) ? region_start_column : "NULL", regions_file);
		}

		if (region_end_column_pos < 0) {
			throw AnnotatorException("Annotator", "process_regions_file_header()", __LINE__, 7, (region_end_column != NULL) ? region_end_column : "NULL", regions_file);
		}
	} catch (ReaderException &e) {
		AnnotatorException new_e(e);
		new_e.add_message("Annotator", "process_regions_file_header()", __LINE__, 6, regions_file);
		throw new_e;
	} catch (DescriptorException &e) {
		AnnotatorException new_e(e);
		new_e.add_message("Annotator", "process_regions_file_header()", __LINE__, 6, regions_file);
		throw new_e;
	}
}

void Annotator::process_regions_file_data()  throw (AnnotatorException) {
	Descriptor* descriptor = NULL;
	char* line = NULL;
	int line_length = 0;
	unsigned int line_number = 2u;

	char data_separator = '\0';

	char* token = NULL;
	char* end_ptr = NULL;

	int column_position = 0;

	const char* region_name_column = NULL;
	const char* region_chr_column = NULL;
	const char* region_start_column = NULL;
	const char* region_end_column = NULL;

	char* region_name_token = NULL;
	char* region_chr_token = NULL;
	char* region_start_token = NULL;
	char* region_end_token = NULL;

	char* region_name = NULL;
	char* region_chr = NULL;
	int region_start = 0;
	int region_end = 0;

	IntervalTree<char*>* regions_index = NULL;

	if ((gwafile == NULL) || (regions_file == NULL)) {
		return;
	}

	try {
		descriptor = gwafile->get_descriptor();
		data_separator = gwafile->get_regions_file_data_separator();

		region_name_column = descriptor->get_property(Descriptor::REGION_NAME);
		region_chr_column = descriptor->get_property(Descriptor::REGION_CHR);
		region_start_column = descriptor->get_property(Descriptor::REGION_START);
		region_end_column = descriptor->get_property(Descriptor::REGION_END);

		while ((line_length = regions_reader.read_line()) > 0) {
			line = *regions_reader.line;

			column_position = 0;
			region_name_token = NULL;
			region_chr_token = NULL;
			region_start_token = NULL;
			region_end_token = NULL;
			token = auxiliary::strtok(&line, data_separator);
			while (token != NULL) {
				if (column_position == region_name_column_pos) {
					region_name_token = token;
				} else if (column_position == region_chr_column_pos) {
					region_chr_token = token;
				} else if (column_position == region_start_column_pos) {
					region_start_token = token;
				} else if (column_position == region_end_column_pos) {
					region_end_token = token;
				}
				token = auxiliary::strtok(&line, data_separator);
				++column_position;
			}

			if (column_position < regions_file_total_columns) {
				throw AnnotatorException("Annotator", "process_regions_file_data()", __LINE__, 8, line_number, regions_file, column_position, regions_file_total_columns);
			} else if (column_position > regions_file_total_columns) {
				throw AnnotatorException("Annotator", "process_regions_file_data()", __LINE__, 9, line_number, regions_file, column_position, regions_file_total_columns);
			}

			region_name = (char*)malloc((strlen(region_name_token) + 1u) * sizeof(char));
			if (region_name == NULL) {
				throw AnnotatorException("Annotator", "process_regions_file_data()",  __LINE__, 2, ((strlen(region_name_token) + 1u) * sizeof(char)));
			}
			strcpy(region_name, region_name_token);

			region_start = (int)strtol(region_start_token, &end_ptr, 10);
			if (*end_ptr != '\0') {
				throw AnnotatorException("Annotator", "process_regions_file_data()",  __LINE__, 10, region_start_token, region_start_column, line_number);
			}

			if (region_start < 0) {
				throw AnnotatorException("Annotator", "process_regions_file_data()",  __LINE__, 11, region_start_column, region_start, line_number);
			}

			region_end = (int)strtol(region_end_token, &end_ptr, 10);
			if (*end_ptr != '\0') {
				throw AnnotatorException("Annotator", "process_regions_file_data()",  __LINE__, 10, region_end_token, region_end_column, line_number);
			}

			if (region_end < 0) {
				throw AnnotatorException("Annotator", "process_regions_file_data()",  __LINE__, 11, region_end_column, region_end, line_number);
			}

			if (region_start > region_end) {
				throw AnnotatorException("Annotator", "process_regions_file_data()",  __LINE__, 12, region_start_column, region_start, region_end, region_end_column, line_number);
			}

			regions_indices_it = regions_indices.find(region_chr_token);
			if (regions_indices_it == regions_indices.end()) {
				regions_index = new IntervalTree<char*>();

				region_chr = (char*)malloc((strlen(region_chr_token) + 1u) * sizeof(char));
				if (region_chr == NULL) {
					throw AnnotatorException("Annotator", "process_regions_file_data()",  __LINE__, 2, ((strlen(region_chr_token) + 1u) * sizeof(char)));
				}
				strcpy(region_chr, region_chr_token);

				regions_indices.insert(pair<char*, IntervalTree<char*>* >(region_chr, regions_index));
			} else {
				regions_index = regions_indices_it->second;
			}

			regions_index->add(region_start, region_end, region_name);
			++line_number;
		}

		if (line_length == 0) {
			throw AnnotatorException("Annotator", "process_regions_file_data()", __LINE__, 13, line_number, regions_file);
		}
	} catch (ReaderException &e) {
		AnnotatorException new_e(e);
		new_e.add_message("Annotator", "process_regions_file_data()", __LINE__, 14, regions_file);
		throw new_e;
	} catch (DescriptorException &e) {
		AnnotatorException new_e(e);
		new_e.add_message("Annotator", "process_regions_file_data()", __LINE__, 14, regions_file);
		throw new_e;
	} catch (AnnotatorException &e) {
		e.add_message("Annotator", "process_regions_file_data()", __LINE__, 14, regions_file);
		throw;
	}
}

void Annotator::write_char_vector(ofstream &ofile_stream, vector<char*>* values, char separator) throw (ofstream::failure) {
	int size = 0;
	if ((values != NULL) && ((size = values->size()) > 0)) {
		--size;
		for (int i = 0; i < size; ++i) {
			ofile_stream << values->at(i) << separator;
		}
		ofile_stream << values->at(size);
	}
}
