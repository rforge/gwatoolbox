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

#ifndef HARMONIZER2_H_
#define HARMONIZER2_H_

#include <cstdio>
#include <map>
#include <set>

#include "Harmonizer2Exception.h"
#include "../../reader/include/ReaderFactory.h"
#include "../../writer/include/WriterFactory.h"
#include "../../auxiliary/include/auxiliary.h"

class Harmonizer2 {
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

	char* vcf_file;
	Reader* vcf_reader;
	unsigned int vcf_file_line_number;
	int vcf_file_column_number;

	char* input_file;
	Reader* input_reader;
	unsigned int input_file_line_number;
	int input_file_column_number;

	char* output_file;
	Writer* output_writer;

	char* log_file;
	Writer* log_writer;

	char** tokens;

	char* id_column;
	char* chr_column;
	char* pos_column;
	char* first_allele_column;
	char* second_allele_column;

	int id_column_pos;
	int chr_column_pos;
	int pos_column_pos;
	int first_allele_column_pos;
	int second_allele_column_pos;

	char separator;
	char* header_backup;

	struct position_index_entry {
		unsigned long int position;
		char* id;
		char type;
		char subtype;
		const char* ref_allele;
		const char* nonref_allele;
	};

	struct chr_index {
		position_index_entry* positions;
		unsigned int n;
		unsigned int heap_n;

		~chr_index() {
			for (unsigned int i = 0u; i < n; ++i) {
				free(positions[i].id);
				positions[i].id = NULL;
				positions[i].ref_allele = NULL;
				positions[i].nonref_allele = NULL;
			}

			free(positions);
			positions = NULL;
		}
	};

	map<char*, chr_index*, bool(*)(const char*, const char*)> index_by_chr;
	map<char*, chr_index*, bool(*)(const char*, const char*)>::iterator index_by_chr_it;

	set<char*, bool(*)(const char*, const char*)> unique_alleles;
	set<char*, bool(*)(const char*, const char*)>::iterator unique_alleles_it;

	static int qsort_position_index_entry_cmp(const void* first, const void* second);

	inline void write_columns() throw (WriterException);

	void open_vcf_file(const char* file_name) throw (Harmonizer2Exception);
	void close_vcf_file() throw (Harmonizer2Exception);
	void process_vcf_file_header() throw (Harmonizer2Exception);
	void process_vcf_file_data() throw (Harmonizer2Exception);

public:
	static const unsigned int INDEX_HEAP_SIZE;
	static const unsigned int INDEX_HEAP_INCREMENT;

	Harmonizer2();
	virtual ~Harmonizer2();

	void open_input_file(const char* file_name, const char* id_column_name, const char* chr_column_name, const char* pos_column_name, const char* first_allele_column_name, const char* second_allele_column_name, char field_separator) throw (Harmonizer2Exception);
	void open_output_file(const char* file_name, bool gzip) throw (Harmonizer2Exception);
	void open_log_file(const char* file_name, bool gzip) throw (Harmonizer2Exception);

	void close_input_file() throw (Harmonizer2Exception);
	void close_output_file() throw (Harmonizer2Exception);
	void close_log_file() throw (Harmonizer2Exception);

	void process_header() throw (Harmonizer2Exception);
	void index_vcf(const char* file_name) throw (Harmonizer2Exception);
	void harmonize(bool flip, bool drop) throw (Harmonizer2Exception);
};

#endif
