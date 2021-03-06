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

#include "include/GwaFile.h"

GwaFile::GwaFile(Descriptor* descriptor) throw (GwaFileException) : descriptor(descriptor) {
	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "GwaFile( Descriptor* )", __LINE__, 0, "descriptor");
	}

	check_filters(descriptor);
	check_thresholds(descriptor);
	check_prefix(descriptor);
	check_casesensitivity(descriptor);
	check_missing_value(descriptor);
	check_separators(descriptor);
	check_filesize(descriptor);
	check_order(descriptor);
	check_genomiccontrol(descriptor);
}

GwaFile::GwaFile(Descriptor* descriptor, void (GwaFile::*check_functions[])(Descriptor*), unsigned int check_functions_number) throw (GwaFileException): descriptor(descriptor)  {
	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "GwaFile( Descriptor* )", __LINE__, 0, "descriptor");
	}

	for (unsigned int i = 0; i < check_functions_number; i++) {
		(this->*check_functions[i])(descriptor);
	}
}

GwaFile::~GwaFile() {
	descriptor = NULL;
}

Descriptor* GwaFile::get_descriptor() {
	return descriptor;
}

bool GwaFile::is_case_sensitive() {
	return casesensitive;
}

bool GwaFile::is_high_verbosity() {
	return high_verbosity;
}

bool GwaFile::is_order_on() {
	return order_on;
}

bool GwaFile::is_gc_on() {
	return gc_on;
}

double GwaFile::get_inflation_factor() {
	return inflation_factor;
}

char GwaFile::get_header_separator() {
	return header_separator;
}

char GwaFile::get_data_separator() {
	return data_separator;
}

unsigned long int GwaFile::get_estimated_size() {
	return estimated_size;
}

char GwaFile::get_regions_file_header_separator() {
	return regions_file_header_separator;
}

char GwaFile::get_regions_file_data_separator() {
	return regions_file_data_separator;
}

bool GwaFile::is_regions_append_on() {
	return regions_append_on;
}

char GwaFile::get_map_file_header_separartor() {
	return map_file_header_separator;
}

char GwaFile::get_map_file_data_separator() {
	return map_file_data_separator;
}

char GwaFile::get_ld_file_header_separator() {
	return ld_file_header_separator;
}

char GwaFile::get_ld_file_data_separator() {
	return ld_file_data_separator;
}

void GwaFile::check_filters(Descriptor* descriptor) throw (GwaFileException) {
	vector<double>* filters = NULL;
	double swap_value = 0.0;

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_filters( Descriptor* )", __LINE__, 0, "descriptor");
	}

	try {
		if ((filters = descriptor->get_threshold(Descriptor::MAF)) == NULL) {
			descriptor->add_threshold(Descriptor::MAF, 0.01);
			descriptor->add_threshold(Descriptor::MAF, 0.05);
		} else if (filters->size() < 2) {
			throw GwaFileException("GwaFile", "check_filters( Descriptor* )", __LINE__, 3, Descriptor::MAF, descriptor->get_full_path());
		} else if (filters->size() > 2) {
			throw GwaFileException("GwaFile", "check_filters( Descriptor* )", __LINE__, 4, Descriptor::MAF, descriptor->get_full_path());
		} else if (filters->at(0) > filters->at(1)) {
			swap_value = filters->at(0);
			filters->at(0) = filters->at(1);
			filters->at(1) = swap_value;
		}

		if ((filters = descriptor->get_threshold(Descriptor::IMP)) == NULL) {
			descriptor->add_threshold(Descriptor::IMP, 0.3);
			descriptor->add_threshold(Descriptor::IMP, 0.6);
		} else if (filters->size() < 2) {
			throw GwaFileException("GwaFile", "check_filters( Descriptor* )", __LINE__, 3, Descriptor::IMP, descriptor->get_full_path());
		} else if (filters->size() > 2) {
			throw GwaFileException("GwaFile", "check_filters( Descriptor* )", __LINE__, 4, Descriptor::IMP, descriptor->get_full_path());
		} else 	if (filters->at(0) > filters->at(1)) {
			swap_value = filters->at(0);
			filters->at(0) = filters->at(1);
			filters->at(1) = swap_value;
		}

		if ((filters = descriptor->get_threshold(Descriptor::SNP_HQ)) == NULL) {
			descriptor->add_threshold(Descriptor::SNP_HQ, 0.01);
			descriptor->add_threshold(Descriptor::SNP_HQ, 0.3);
		} else if (filters->size() < 2) {
			throw GwaFileException("GwaFile", "check_filters( Descriptor* )", __LINE__, 5, Descriptor::SNP_HQ, descriptor->get_full_path());
		} else if (filters->size() > 2) {
			throw GwaFileException("GwaFile", "check_filters( Descriptor* )", __LINE__, 6, Descriptor::SNP_HQ, descriptor->get_full_path());
		}
	} catch (DescriptorException &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_filters( Descriptor* )", __LINE__, 10, descriptor->get_full_path());
		throw new_e;
	}
}

void GwaFile::check_re_filters(Descriptor* descriptor) throw (GwaFileException) {
	vector<double>* filters = NULL;
	double swap_value = 0.0;

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_filters( Descriptor* )", __LINE__, 0, "descriptor");
	}

	try {
		if ((filters = descriptor->get_threshold(Descriptor::MAF)) == NULL) {
			descriptor->add_threshold(Descriptor::MAF, 0.01);
			descriptor->add_threshold(Descriptor::MAF, 0.05);
		} else if (filters->size() < 2) {
			throw GwaFileException("GwaFile", "check_re_filters( Descriptor* )", __LINE__, 3, Descriptor::MAF, descriptor->get_full_path());
		} else if (filters->size() > 2) {
			throw GwaFileException("GwaFile", "check_re_filters( Descriptor* )", __LINE__, 4, Descriptor::MAF, descriptor->get_full_path());
		} else if (filters->at(0) > filters->at(1)) {
			swap_value = filters->at(0);
			filters->at(0) = filters->at(1);
			filters->at(1) = swap_value;
		}

		if ((filters = descriptor->get_threshold(Descriptor::IMP)) == NULL) {
			descriptor->add_threshold(Descriptor::IMP, 0.3);
		} else if (filters->size() > 1) {
			throw GwaFileException("GwaFile", "check_re_filters( Descriptor* )", __LINE__, 16, Descriptor::IMP, descriptor->get_full_path());
		}
	} catch (DescriptorException &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_randomeffect_filters( Descriptor* )", __LINE__, 10, descriptor->get_full_path());
		throw new_e;
	}
}

void GwaFile::check_thresholds(Descriptor* descriptor) throw (GwaFileException) {
	vector<double>* thresholds = NULL;
	double swap_value = 0.0;

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 0, "descriptor");
	}

	try {
		if ((thresholds = descriptor->get_threshold(Descriptor::STDERR)) == NULL) {
			descriptor->add_threshold(Descriptor::STDERR, 0.0);
			descriptor->add_threshold(Descriptor::STDERR, 100000.0);
		} else if (thresholds->size() < 2) {
			throw GwaFileException("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 7, Descriptor::STDERR, descriptor->get_full_path(), 2);
		} else if (thresholds->size() > 2) {
			throw GwaFileException("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 8, Descriptor::STDERR, descriptor->get_full_path(), 2);
		} else if (thresholds->at(0) > thresholds->at(1)) {
			swap_value = thresholds->at(0);
			thresholds->at(0) = thresholds->at(1);
			thresholds->at(1) = swap_value;
		}

		if ((thresholds = descriptor->get_threshold(Descriptor::OEVAR_IMP)) == NULL) {
			descriptor->add_threshold(Descriptor::OEVAR_IMP, 0.0);
			descriptor->add_threshold(Descriptor::OEVAR_IMP, 1.5);
		} else if (thresholds->size() < 2) {
			throw GwaFileException("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 7, Descriptor::OEVAR_IMP, descriptor->get_full_path(), 2);
		} else if (thresholds->size() > 2) {
			throw GwaFileException("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 8, Descriptor::OEVAR_IMP, descriptor->get_full_path(), 2);
		} else if (thresholds->at(0) > thresholds->at(1)) {
			swap_value = thresholds->at(0);
			thresholds->at(0) = thresholds->at(1);
			thresholds->at(1) = swap_value;
		}

		if ((thresholds = descriptor->get_threshold(Descriptor::PVALUE)) == NULL) {
			descriptor->add_threshold(Descriptor::PVALUE, 0.0);
			descriptor->add_threshold(Descriptor::PVALUE, 1.0);
		} else if (thresholds->size() < 2) {
			throw GwaFileException("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 7, Descriptor::PVALUE, descriptor->get_full_path(), 2);
		} else if (thresholds->size() > 2) {
			throw GwaFileException("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 8, Descriptor::PVALUE, descriptor->get_full_path(), 2);
		} else if (thresholds->at(0) > thresholds->at(1)) {
			swap_value = thresholds->at(0);
			thresholds->at(0) = thresholds->at(1);
			thresholds->at(1) = swap_value;
		}

		if ((thresholds = descriptor->get_threshold(Descriptor::FREQLABEL)) == NULL) {
			descriptor->add_threshold(Descriptor::FREQLABEL, 0.0);
			descriptor->add_threshold(Descriptor::FREQLABEL, 1.0);
		} else if (thresholds->size() < 2) {
			throw GwaFileException("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 7, Descriptor::FREQLABEL, descriptor->get_full_path(), 2);
		} else if (thresholds->size() > 2) {
			throw GwaFileException("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 8, Descriptor::FREQLABEL, descriptor->get_full_path(), 2);
		} else if (thresholds->at(0) > thresholds->at(1)) {
			swap_value = thresholds->at(0);
			thresholds->at(0) = thresholds->at(1);
			thresholds->at(1) = swap_value;
		}

		if ((thresholds = descriptor->get_threshold(Descriptor::HWE_PVAL)) == NULL) {
			descriptor->add_threshold(Descriptor::HWE_PVAL, 0.0);
			descriptor->add_threshold(Descriptor::HWE_PVAL, 1.0);
		} else if (thresholds->size() < 2) {
			throw GwaFileException("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 7, Descriptor::HWE_PVAL, descriptor->get_full_path(), 2);
		} else if (thresholds->size() > 2) {
			throw GwaFileException("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 8, Descriptor::HWE_PVAL, descriptor->get_full_path(), 2);
		} else if (thresholds->at(0) > thresholds->at(1)) {
			swap_value = thresholds->at(0);
			thresholds->at(0) = thresholds->at(1);
			thresholds->at(1) = swap_value;
		}

		if ((thresholds = descriptor->get_threshold(Descriptor::CALLRATE)) == NULL) {
			descriptor->add_threshold(Descriptor::CALLRATE, 0.0);
			descriptor->add_threshold(Descriptor::CALLRATE, 1.0);
		} else if (thresholds->size() < 2) {
			throw GwaFileException("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 7, Descriptor::CALLRATE, descriptor->get_full_path(), 2);
		} else if (thresholds->size() > 2) {
			throw GwaFileException("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 8, Descriptor::CALLRATE, descriptor->get_full_path(), 2);
		} else if (thresholds->at(0) > thresholds->at(1)) {
			swap_value = thresholds->at(0);
			thresholds->at(0) = thresholds->at(1);
			thresholds->at(1) = swap_value;
		}
	} catch (DescriptorException &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 10, descriptor->get_full_path());
		throw new_e;
	}
}

void GwaFile::check_prefix(Descriptor* descriptor) throw (GwaFileException) {
	const char* prefix = NULL;

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_prefix( Descriptor* )", __LINE__, 0, "descriptor");
	}

	try {
		if ((prefix = descriptor->get_property(Descriptor::PREFIX)) == NULL) {
			throw GwaFileException("GwaFile", "check_prefix( Descriptor* )", __LINE__, 9, Descriptor::PREFIX, descriptor->get_full_path());
		}
	} catch (DescriptorException &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_prefix( Descriptor* )", __LINE__, 10, descriptor->get_full_path());
		throw new_e;
	}
}

void GwaFile::check_casesensitivity(Descriptor* descriptor) throw (GwaFileException) {
	const char* casesensitive = NULL;

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_casesensitivity( Descriptor* )", __LINE__, 0, "descriptor");
	}

	try {
		this->casesensitive = false;
		if ((casesensitive = descriptor->get_property(Descriptor::CASESENSITIVE)) == NULL) {
			descriptor->add_property(Descriptor::CASESENSITIVE, Descriptor::OFF_MODES[0]);
		} else if ((strcmp_ignore_case(casesensitive, Descriptor::ON_MODES[0]) == 0) ||
				(strcmp_ignore_case(casesensitive, Descriptor::ON_MODES[1]) == 0)) {
			this->casesensitive = true;
		}
	} catch (DescriptorException &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_casesensitivity( Descriptor* )", __LINE__, 10, descriptor->get_full_path());
		throw new_e;
	}
}

void GwaFile::check_verbosity_level(Descriptor* descriptor) throw (GwaFileException) {
	const char* verbosity_level = NULL;

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_verbosity_level( Descriptor* )", __LINE__, 0, "descriptor");
	}

	try {
		high_verbosity = false;
		if ((verbosity_level = descriptor->get_property(Descriptor::VERBOSITY)) == NULL) {
			descriptor->add_property(Descriptor::VERBOSITY, Descriptor::LOW_MODES[0]);
		} else if ((strcmp_ignore_case(verbosity_level, Descriptor::HIGH_MODES[0]) == 0) ||
				(strcmp_ignore_case(verbosity_level, Descriptor::HIGH_MODES[1]) == 0)) {
			high_verbosity = true;
		}
	} catch (DescriptorException &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_verbosity_level( Descriptor* )", __LINE__, 10, descriptor->get_full_path());
		throw new_e;
	}
}

void GwaFile::check_missing_value(Descriptor* descriptor) throw (GwaFileException) {
	const char* missing_value = NULL;

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_missing_value( Descriptor* )", __LINE__, 0, "descriptor");
	}

	try {
		if ((missing_value = descriptor->get_property(Descriptor::MISSING)) == NULL) {
			descriptor->add_property(Descriptor::MISSING, "NA");
		}
	} catch (DescriptorException &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_missing_value( Descriptor* )", __LINE__, 10, descriptor->get_full_path());
		throw new_e;
	}
}

void GwaFile::check_separators(Descriptor* descriptor) throw (GwaFileException) {
	const char* separator_name = NULL;

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_separators( Descriptor* )", __LINE__, 0, "descriptor");
	}

	try {
		if ((separator_name = descriptor->get_property(Descriptor::SEPARATOR)) == NULL) {
			Reader* reader = NULL;

			reader = ReaderFactory::create(descriptor->get_full_path());
			reader->detect_field_separators(&header_separator, &data_separator);
			reader->close();

			delete reader;
			reader = NULL;
		} else if ((strcmp_ignore_case(separator_name, Descriptor::COMMA) == 0) ||
				(strcmp_ignore_case(separator_name, Descriptor::COMMAS) == 0)) {
			header_separator = ',';
			data_separator = ',';
		} else if (strcmp_ignore_case(separator_name, Descriptor::SEMICOLON) == 0) {
			header_separator = ';';
			data_separator = ';';
		} else if ((strcmp_ignore_case(separator_name, Descriptor::TAB) == 0) ||
				(strcmp_ignore_case(separator_name, Descriptor::TABULATION) == 0)) {
			header_separator = '\t';
			data_separator = '\t';
		} else if (strcmp_ignore_case(separator_name, Descriptor::WHITESPACE) == 0) {
			header_separator = ' ';
			data_separator = ' ';
		}
	} catch (Exception &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_separators( Descriptor* )", __LINE__, 10, descriptor->get_full_path());
		throw new_e;
	}
}

void GwaFile::check_filesize(Descriptor* descriptor) throw (GwaFileException) {
	try {
		TextReader reader;
		reader.set_file_name(descriptor->get_full_path());
		estimated_size = reader.estimate_lines_count();
		reader.close();
	} catch (Exception &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_filesize( Descriptor* )", __LINE__, 11, descriptor->get_full_path());
		throw new_e;
	}
}

void GwaFile::check_order(Descriptor* descriptor) throw (GwaFileException) {
	const char* order_on = NULL;

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_order( Descriptor* )", __LINE__, 0, "descriptor");
	}

	try {
		this->order_on = false;
		if ((order_on = descriptor->get_property(Descriptor::ORDER)) == NULL) {
			descriptor->add_property(Descriptor::ORDER, Descriptor::OFF_MODES[0]);
		} else if ((strcmp_ignore_case(order_on, Descriptor::ON_MODES[0]) == 0) ||
				(strcmp_ignore_case(order_on, Descriptor::ON_MODES[1]) == 0)) {
			this->order_on = true;
		}
	} catch (DescriptorException &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_order( Descriptor* )", __LINE__, 10, descriptor->get_full_path());
		throw new_e;
	}
}

void GwaFile::check_genomiccontrol(Descriptor* descriptor) throw (GwaFileException) {
	const char* gc = NULL;

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_genomiccontrol( Descriptor* )", __LINE__, 0, "descriptor");
	}

	try {
		inflation_factor = numeric_limits<double>::quiet_NaN();
		if ((gc = descriptor->get_property(Descriptor::GENOMICCONTROL)) == NULL) {
			descriptor->add_property(Descriptor::GENOMICCONTROL, Descriptor::OFF_MODES[0]);
			gc_on = false;
		} else if (strcmp_ignore_case(gc, Descriptor::ON_MODES[0]) == 0) {
			gc_on = true;
		} else if (strcmp_ignore_case(gc, Descriptor::OFF_MODES[0]) == 0) {
			gc_on = false;
		} else {
			inflation_factor = auxiliary::to_double(gc);
			gc_on = true;
		}
	} catch (DescriptorException &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_order( Descriptor* )", __LINE__, 10, descriptor->get_full_path());
		throw new_e;
	}
}

void GwaFile::check_regions_file(Descriptor* descriptor) throw (GwaFileException) {
	const char* genes_file = NULL;

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_regions_file( Descriptor* )", __LINE__, 0, "descriptor");
	}

	try {
		if ((genes_file = descriptor->get_property(Descriptor::REGIONS_FILE)) == NULL) {
			throw GwaFileException("GwaFile", "check_regions_file( Descriptor* )", __LINE__, 12, Descriptor::REGIONS_FILE, descriptor->get_full_path());
		}
	} catch (DescriptorException &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_regions_file( Descriptor* )", __LINE__, 10, descriptor->get_full_path());
		throw new_e;
	}
}

void GwaFile::check_regions_file_separators(Descriptor* descriptor) throw (GwaFileException) {
	const char* genes_file = NULL;
	const char* separator_name = NULL;

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_regions_file_separators( Descriptor* )", __LINE__, 0, "descriptor");
	}

	try {
		if ((genes_file = descriptor->get_property(Descriptor::REGIONS_FILE)) == NULL) {
			throw GwaFileException("GwaFile", "check_regions_file_separators( Descriptor* )", __LINE__, 12, Descriptor::REGIONS_FILE, descriptor->get_full_path());
		}

		if ((separator_name = descriptor->get_property(Descriptor::REGIONS_FILE_SEPARATOR)) == NULL) {
			TextReader reader;
			reader.set_file_name(genes_file);
			reader.detect_field_separators(&regions_file_header_separator, &regions_file_data_separator);
			reader.close();
		} else if ((strcmp_ignore_case(separator_name, Descriptor::COMMA) == 0) ||
				(strcmp_ignore_case(separator_name, Descriptor::COMMAS) == 0)) {
			regions_file_header_separator = ',';
			regions_file_data_separator = ',';
		} else if (strcmp_ignore_case(separator_name, Descriptor::SEMICOLON) == 0) {
			regions_file_header_separator = ';';
			regions_file_data_separator = ';';
		} else if ((strcmp_ignore_case(separator_name, Descriptor::TAB) == 0) ||
				(strcmp_ignore_case(separator_name, Descriptor::TABULATION) == 0)) {
			regions_file_header_separator = '\t';
			regions_file_data_separator = '\t';
		} else if (strcmp_ignore_case(separator_name, Descriptor::WHITESPACE) == 0) {
			regions_file_header_separator = ' ';
			regions_file_data_separator = ' ';
		}
	} catch (Exception &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_regions_file_separators( Descriptor* )", __LINE__, 10, descriptor->get_full_path());
		throw new_e;
	}
}

void GwaFile::check_regions_deviation(Descriptor* descriptor) throw (GwaFileException) {
	vector<double>* deviation = NULL;

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_regions_deviation( Descriptor* )", __LINE__, 0, "descriptor");
	}

	if ((deviation = descriptor->get_threshold(Descriptor::REGIONS_DEVIATION)) == NULL) {
		descriptor->add_threshold(Descriptor::REGIONS_DEVIATION, 0.0);
		descriptor->add_threshold(Descriptor::REGIONS_DEVIATION, 100000.0);
		descriptor->add_threshold(Descriptor::REGIONS_DEVIATION, 250000.0);
	} else {
		sort(deviation->begin(), deviation->end());
		if (deviation->at(0) < 0.0) {
			throw GwaFileException("GwaFile", "check_regions_deviation( Descriptor* )", __LINE__, 13, Descriptor::REGIONS_DEVIATION, descriptor->get_full_path());
		}
	}
}

void GwaFile::check_regions_append(Descriptor* descriptor) throw (GwaFileException) {
	const char* regions_append = NULL;

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_regions_append( Descriptor* )", __LINE__, 0, "descriptor");
	}

	try {
		regions_append_on = false;
		if ((regions_append = descriptor->get_property(Descriptor::REGIONS_APPEND)) == NULL) {
			descriptor->add_property(Descriptor::REGIONS_APPEND, Descriptor::OFF_MODES[0]);
		} else if ((strcmp_ignore_case(regions_append, Descriptor::ON_MODES[0]) == 0) ||
				(strcmp_ignore_case(regions_append, Descriptor::ON_MODES[1]) == 0)) {
			regions_append_on = true;
		}
	} catch (DescriptorException &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_regions_append( Descriptor* )", __LINE__, 10, descriptor->get_full_path());
		throw new_e;
	}
}

void GwaFile::check_map_file_separators(Descriptor* descriptor) throw (GwaFileException) {
	const char* map_file = NULL;
	const char* separator_name = NULL;

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_map_file_separators( Descriptor* )", __LINE__, 0, "descriptor");
	}

	try {
		if ((map_file = descriptor->get_property(Descriptor::MAP_FILE)) != NULL) {
			if ((separator_name = descriptor->get_property(Descriptor::MAP_FILE_SEPARATOR)) == NULL) {
				TextReader reader;
				reader.set_file_name(map_file);
				reader.detect_field_separators(&map_file_header_separator, &map_file_data_separator);
				reader.close();
			} else if ((strcmp_ignore_case(separator_name, Descriptor::COMMA) == 0) ||
					(strcmp_ignore_case(separator_name, Descriptor::COMMAS) == 0)) {
				map_file_header_separator = ',';
				map_file_data_separator = ',';
			} else if (strcmp_ignore_case(separator_name, Descriptor::SEMICOLON) == 0) {
				map_file_header_separator = ';';
				map_file_data_separator = ';';
			} else if ((strcmp_ignore_case(separator_name, Descriptor::TAB) == 0) ||
					(strcmp_ignore_case(separator_name, Descriptor::TABULATION) == 0)) {
				map_file_header_separator = '\t';
				map_file_data_separator = '\t';
			} else if (strcmp_ignore_case(separator_name, Descriptor::WHITESPACE) == 0) {
				map_file_header_separator = ' ';
				map_file_data_separator = ' ';
			}
		}
	} catch (Exception &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_map_file_separators( Descriptor* )", __LINE__, 10, descriptor->get_full_path());
		throw new_e;
	}
}

void GwaFile::check_ld_files(Descriptor* descriptor) throw (GwaFileException) {
	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_ld_files( Descriptor* )", __LINE__, 0, "descriptor");
	}

	if (descriptor->get_ld_files_number() <= 0) {
		throw GwaFileException("GwaFile", "check_ld_files( Descriptor )", __LINE__, 14, Descriptor::LD_FILE, descriptor->get_full_path());
	}
}

void GwaFile::check_ld_files_separators(Descriptor* descriptor) throw (GwaFileException) {
	const char* separator_name = NULL;
	vector<const char*>* files = NULL;
	vector<const char*>::iterator files_it;
	Reader* reader = NULL;

	char ld_file_header_separator = '\0';
	char ld_file_data_separator = '\0';

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_ld_files_separators( Descriptor* )", __LINE__, 0, "descriptor");
	}

	try {
		if ((separator_name = descriptor->get_property(Descriptor::LD_FILE_SEPARATOR)) == NULL) {
			files = descriptor->get_ld_files();
			if ((files == NULL) || (files->size() <= 0)) {
				throw GwaFileException("GwaFile", "check_ld_files_separators( Descriptor* )", __LINE__, 14, Descriptor::LD_FILE, descriptor->get_full_path());
			}

			files_it = files->begin();

			reader = ReaderFactory::create(*files_it);
			reader->detect_field_separators(&(this->ld_file_header_separator), &(this->ld_file_data_separator));
			reader->close();
			delete reader;
			reader = NULL;

			files_it++;

			while (files_it != files->end()) {
				reader = ReaderFactory::create(*files_it);
				reader->detect_field_separators(&ld_file_header_separator, &ld_file_data_separator);
				reader->close();
				delete reader;
				reader = NULL;

				if ((ld_file_header_separator != this->ld_file_header_separator) || (ld_file_data_separator != this->ld_file_data_separator)) {
					throw GwaFileException("GwaFile", "check_ld_files_separators( Descriptor* )", __LINE__, 15, *files_it, Descriptor::LD_FILE, descriptor->get_full_path());
				}

				files_it++;
			}

			delete files;
			files = NULL;
		} else if ((strcmp_ignore_case(separator_name, Descriptor::COMMA) == 0) ||
				(strcmp_ignore_case(separator_name, Descriptor::COMMAS) == 0)) {
			this->ld_file_header_separator = ',';
			this->ld_file_data_separator = ',';
		} else if (strcmp_ignore_case(separator_name, Descriptor::SEMICOLON) == 0) {
			this->ld_file_header_separator = ';';
			this->ld_file_data_separator = ';';
		} else if ((strcmp_ignore_case(separator_name, Descriptor::TAB) == 0) ||
				(strcmp_ignore_case(separator_name, Descriptor::TABULATION) == 0)) {
			this->ld_file_header_separator = '\t';
			this->ld_file_data_separator = '\t';
		} else if (strcmp_ignore_case(separator_name, Descriptor::WHITESPACE) == 0) {
			this->ld_file_header_separator = ' ';
			this->ld_file_data_separator = ' ';
		}
	} catch (Exception &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_ld_files_separators( Descriptor* )", __LINE__, 10, descriptor->get_full_path());
		throw new_e;
	}
}

void GwaFile::check_ld_threshold(Descriptor* descriptor) throw (GwaFileException) {
	vector<double>* thresholds = NULL;

	if (descriptor == NULL) {
		throw GwaFileException("GwaFile", "check_ld_threshold( Descriptor* )", __LINE__, 0, "descriptor");
	}

	try {
		if ((thresholds = descriptor->get_threshold(Descriptor::LD)) == NULL) {
			descriptor->add_threshold(Descriptor::LD, 0.3);
		} else if (thresholds->size() < 1) {
			throw GwaFileException("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 7, Descriptor::LD, descriptor->get_full_path(), 1);
		} else if (thresholds->size() > 1) {
			throw GwaFileException("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 8, Descriptor::LD, descriptor->get_full_path(), 1);
		}
	} catch (DescriptorException &e) {
		GwaFileException new_e(e);
		new_e.add_message("GwaFile", "check_thresholds( Descriptor* )", __LINE__, 10, descriptor->get_full_path());
		throw new_e;
	}
}
