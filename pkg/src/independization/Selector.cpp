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

#include "include/Selector.h"

Selector::Selector() : gwafile(NULL), reader(NULL), header_backup(NULL),
	total_columns(numeric_limits<int>::min()),
	marker_column_pos(numeric_limits<int>::min()),
	chr_column_pos(numeric_limits<int>::min()),
	pvalue_column_pos(numeric_limits<int>::min()) {

}

Selector::~Selector() {
	gwafile = NULL;

	if (reader != NULL) {
		reader->close();

		delete reader;
		reader = NULL;
	}

	if (header_backup != NULL) {
		free(header_backup);
		header_backup = NULL;
	}
}

void Selector::open_gwafile(GwaFile* gwafile) throw (SelectorException) {
	if (gwafile == NULL) {
		throw SelectorException("Selector", "open_gwafile( GwaFile* )", __LINE__, 0, "gwafile");
	}

	try {
		close_gwafile();

		this->gwafile = gwafile;

		reader = ReaderFactory::create(gwafile->get_descriptor()->get_full_path());
		reader->open();
	} catch (DescriptorException& e) {
		SelectorException new_e(e);
		new_e.add_message("Selector", "open_gwafile( GwaFile* )", __LINE__, 3, gwafile->get_descriptor()->get_full_path());
		throw new_e;
	} catch (ReaderException& e) {
		SelectorException new_e(e);
		new_e.add_message("Selector", "open_gwafile( GwaFile* )", __LINE__, 3, gwafile->get_descriptor()->get_full_path());
		throw new_e;
	} catch (SelectorException& e) {
		e.add_message("Selector", "open_gwafile( GwaFile* )", __LINE__, 3, gwafile->get_descriptor()->get_full_path());
		throw;
	}
}

void Selector::close_gwafile() throw (SelectorException) {
	try {
		if (reader != NULL) {
			reader->close();
			delete reader;
			reader = NULL;
		}
	} catch (ReaderException &e) {
		SelectorException new_e(e);
		new_e.add_message("Selector", "close_gwafile()", __LINE__, 4, (gwafile != NULL) ? gwafile->get_descriptor()->get_full_path() : "NULL");
		throw new_e;
	}

	if (header_backup != NULL) {
		free(header_backup);
		header_backup = NULL;
	}

	gwafile = NULL;
}

void Selector::process_header() throw (SelectorException) {
	Descriptor* descriptor = NULL;
	char header_separator = '\0';
	char* header = NULL;
	char* token = NULL;
	int column_position = 0;
	const char* column_name = NULL;

	if (gwafile == NULL) {
		return;
	}

	try {
		descriptor = gwafile->get_descriptor();
		header_separator = gwafile->get_header_separator();

		if (reader->read_line() <= 0) {
			throw SelectorException("Selector", "process_header()", __LINE__, 5, 1, gwafile->get_descriptor()->get_name());
		}

		header = *(reader->line);

		header_backup = (char*)malloc((strlen(header) + 1u) * sizeof(char));
		if (header_backup == NULL) {
			throw SelectorException("Selector", "process_header()", __LINE__, 2, ((strlen(header) + 1u) * sizeof(char)));
		}
		strcpy(header_backup, header);

		total_columns = numeric_limits<int>::min();
		marker_column_pos = numeric_limits<int>::min();
		chr_column_pos = numeric_limits<int>::min();
		pvalue_column_pos = numeric_limits<int>::min();

		token = auxiliary::strtok(&header, header_separator);
		while (token != NULL) {
			column_name = descriptor->get_default_column(token, gwafile->is_case_sensitive());
			if (column_name != NULL) {
				if (strcmp(column_name, Descriptor::MARKER) == 0) {
					marker_column_pos = column_position;
				} else if (strcmp(column_name, Descriptor::CHR) == 0) {
					chr_column_pos = column_position;
				} else if (strcmp(column_name, Descriptor::PVALUE) == 0) {
					pvalue_column_pos = column_position;
				}
			}
			token = auxiliary::strtok(&header, header_separator);
			++column_position;
		}

		total_columns = column_position;

		if (marker_column_pos < 0) {
			throw SelectorException("Selector", "process_header()", __LINE__, 7, ((column_name = descriptor->get_column(Descriptor::MARKER)) != NULL) ? column_name : Descriptor::MARKER, gwafile->get_descriptor()->get_name());
		}

		if (chr_column_pos < 0) {
			throw SelectorException("Selector", "process_header()", __LINE__, 7, ((column_name = descriptor->get_column(Descriptor::CHR)) != NULL) ? column_name : Descriptor::CHR, gwafile->get_descriptor()->get_name());
		}

		if (pvalue_column_pos < 0) {
			throw SelectorException("Selector", "process_header()", __LINE__, 7, ((column_name = descriptor->get_column(Descriptor::PVALUE)) != NULL) ? column_name : Descriptor::PVALUE, gwafile->get_descriptor()->get_name());
		}
	} catch (ReaderException &e) {
		SelectorException new_e(e);
		new_e.add_message("Selector", "process_header()", __LINE__, 6, gwafile->get_descriptor()->get_name());
		throw new_e;
	} catch (DescriptorException &e) {
		SelectorException new_e(e);
		new_e.add_message("Selector", "process_header()", __LINE__, 6, gwafile->get_descriptor()->get_name());
		throw new_e;
	}
}
