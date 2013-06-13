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

#ifndef HARMONIZER_H_
#define HARMONIZER_H_

#include <cstdio>
#include <map>

#include "HarmonizerException.h"
#include "../../reader/include/ReaderFactory.h"
#include "../../writer/include/WriterFactory.h"
#include "../../auxiliary/include/auxiliary.h"

class Harmonizer {
private:
	static const char VCF_FIELD_SEPARATOR;
	static const char* VCF_FILE_FORMAT;
	static const char* VCF_CHROM;
	static const char* VCF_POS;
	static const char* VCF_ID;
	static const char* VCF_REF;
	static const char* VCF_ALT;
	static const char* VCF_QUAL;
	static const char* VCF_FILTER;
	static const char* VCF_INFO;
	static const char* VCF_FORMAT;
	static const int VCF_MANDATORY_COLUMNS_SIZE;
	static const char* vcf_mandatory_columns[];
	static const char* VCF_PASS;
	static const char* VCF_MISSING;
	static const char VCF_INFO_FIELD_SEPARATOR;
	static const char* VCF_VARIANT_TYPE;
	static const char* VCF_SNP_TYPE;
	static const char* VCF_INDEL_TYPE_01;
	static const char* VCF_INDEL_TYPE_02;
	static const char* VCF_ALT_ALLELE_DEL;
	static const char* VCF_ALT_ALLELE_INS;

	char* map_file;
	Reader* map_reader;
	unsigned int map_file_line_number;
	int map_file_column_number;

	char* input_file;
	char* output_file;
	char* log_file;

	Reader* reader;
	Writer* writer;
	Writer* log_writer;

	char** tokens;

	char* chr_column;
	char* id_column;
	char* ref_allele_column;
	char* nonref_allele_column;
	char separator;
	char* header_backup;
	int file_column_number;
	int chr_column_pos;
	int id_column_pos;
	int ref_allele_column_pos;
	int nonref_allele_column_pos;

	struct position_index_entry {
		unsigned long int position;
		char* id;
		char type;
		char* ref_allele;
		char* nonref_allele;
	};

	struct id_index_entry {
		const char* id;
		char type;
		unsigned int location;
	};

	struct chr_index {
		position_index_entry* positions;
		id_index_entry* ids;
		unsigned int n;
		unsigned int heap_n;

		~chr_index() {
			for (unsigned int i = 0u; i < n; ++i) {
				free(positions[i].id);
				positions[i].id = NULL;
				if (positions[i].ref_allele != NULL) {
					free(positions[i].ref_allele);
					positions[i].ref_allele = NULL;
				}
				if (positions[i].nonref_allele != NULL) {
					free(positions[i].nonref_allele);
					positions[i].nonref_allele = NULL;
				}
				ids[i].id = NULL;
			}

			free(positions);
			free(ids);

			positions = NULL;
			ids = NULL;
		}
	};

	map<char*, chr_index*, bool(*)(const char*, const char*)> map_index_by_chr;
	map<char*, chr_index*, bool(*)(const char*, const char*)>::iterator map_index_by_chr_it;

	static int qsort_position_index_entry_cmp(const void* first, const void* second);
	static int qsort_id_index_entry_cmp(const void* first, const void* second);

	inline void write_columns() throw (WriterException);

	void open_map_file(const char* file_name) throw (HarmonizerException);
	void close_map_file() throw (HarmonizerException);
	void process_map_file_header() throw (HarmonizerException);
	void process_map_file_data() throw (HarmonizerException);

public:
	static const unsigned int MAP_HEAP_SIZE;
	static const unsigned int MAP_HEAP_INCREMENT;

	Harmonizer();
	virtual ~Harmonizer();

	void open_input_file(const char* file_name, const char* chr_column_name, const char* id_column_name, const char* ref_allele_column_name, const char* nonref_allele_column_name, char field_separator) throw (HarmonizerException);
	void open_output_file(const char* file_name, bool gzip) throw (HarmonizerException);
	void open_log_file(const char* file_name, bool gzip) throw (HarmonizerException);

	void close_input_file() throw (HarmonizerException);
	void close_output_file() throw (HarmonizerException);
	void close_log_file() throw (HarmonizerException);

	void process_header() throw (HarmonizerException);
	void index_map(const char* file_name) throw (HarmonizerException);
	void harmonize(bool drop) throw (HarmonizerException);
	void harmonize_no_vcf_allele_check(bool drop) throw (HarmonizerException);
};

#endif
