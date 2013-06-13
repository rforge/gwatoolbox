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

#include "include/Harmonizer2.h"

const char Harmonizer2::VCF_FIELD_SEPARATOR = '\t';
const char* Harmonizer2::VCF_FILE_FORMAT = "##fileformat";
const char* Harmonizer2::VCF_CHROM = "#CHROM";
const char* Harmonizer2::VCF_POS = "POS";
const char* Harmonizer2::VCF_ID = "ID";
const char* Harmonizer2::VCF_REF = "REF";
const char* Harmonizer2::VCF_ALT = "ALT";
const char* Harmonizer2::VCF_QUAL = "QUAL";
const char* Harmonizer2::VCF_FILTER = "FILTER";
const char* Harmonizer2::VCF_INFO = "INFO";
const char* Harmonizer2::VCF_FORMAT = "FORMAT";
const int Harmonizer2::VCF_MANDATORY_COLUMNS_SIZE = 9;
const char* Harmonizer2::vcf_mandatory_columns[VCF_MANDATORY_COLUMNS_SIZE] = {
		VCF_CHROM, VCF_POS, VCF_ID, VCF_REF, VCF_ALT, VCF_QUAL, VCF_FILTER, VCF_INFO, VCF_FORMAT
};
const char* Harmonizer2::VCF_PASS = "PASS";
const char* Harmonizer2::VCF_MISSING = ".";
const char Harmonizer2::VCF_INFO_FIELD_SEPARATOR = ';';
const char* Harmonizer2::VCF_VARIANT_TYPE = "VT";
const char* Harmonizer2::VCF_SNP_TYPE = "SNP";
const char* Harmonizer2::VCF_INDEL_TYPE_01 = "INDEL";
const char* Harmonizer2::VCF_INDEL_TYPE_02 = "I";
const char* Harmonizer2::VCF_ALT_ALLELE_DEL = "<DEL>";

const unsigned int Harmonizer2::INDEX_HEAP_SIZE = 3000000;
const unsigned int Harmonizer2::INDEX_HEAP_INCREMENT = 1000000;

Harmonizer2::Harmonizer2() :
		vcf_file(NULL), vcf_reader(NULL), vcf_file_line_number(0u), vcf_file_column_number(0),
		input_file(NULL), input_reader(NULL), input_file_line_number(0u), input_file_column_number(0),
		output_file(NULL), output_writer(NULL), log_file(NULL), log_writer(NULL),
		tokens(NULL), id_column(NULL), chr_column(NULL), pos_column(NULL), first_allele_column(NULL), second_allele_column(NULL),
		id_column_pos(numeric_limits<int>::min()), chr_column_pos(numeric_limits<int>::min()), pos_column_pos(numeric_limits<int>::min()),
		first_allele_column_pos(numeric_limits<int>::min()), second_allele_column_pos(numeric_limits<int>::min()),
		separator('\0'), header_backup(NULL),
		index_by_chr(auxiliary::bool_strcmp_ignore_case), unique_alleles(auxiliary::bool_strcmp_ignore_case) {
}

Harmonizer2::~Harmonizer2() {
	if (vcf_file != NULL) {
		free(vcf_file);
		vcf_file = NULL;
	}

	if (vcf_reader != NULL) {
		delete vcf_reader;
		vcf_reader = NULL;
	}

	if (input_file != NULL) {
		free(input_file);
		input_file = NULL;
	}

	if (input_reader != NULL) {
		delete input_reader;
		input_reader = NULL;
	}

	if (output_file != NULL) {
		free(output_file);
		output_file = NULL;
	}

	if (output_writer != NULL) {
		delete output_writer;
		output_writer = NULL;
	}

	if (log_file != NULL) {
		free(log_file);
		log_file = NULL;
	}

	if (log_writer != NULL) {
		delete log_writer;
		log_writer = NULL;
	}

	if (tokens != NULL) {
		delete tokens;
		tokens = NULL;
	}

	if (id_column != NULL) {
		free(id_column);
		id_column = NULL;
	}

	if (chr_column != NULL) {
		free(chr_column);
		chr_column = NULL;
	}

	if (pos_column != NULL) {
		free(pos_column);
		pos_column = NULL;
	}

	if (first_allele_column != NULL) {
		free(first_allele_column);
		first_allele_column = NULL;
	}

	if (second_allele_column != NULL) {
		free(second_allele_column);
		second_allele_column = NULL;
	}

	if (header_backup != NULL) {
		free(header_backup);
		header_backup = NULL;
	}

	for (index_by_chr_it = index_by_chr.begin(); index_by_chr_it != index_by_chr.end(); ++index_by_chr_it) {
		free(index_by_chr_it->first);
		delete index_by_chr_it->second;
	}

	for (unique_alleles_it = unique_alleles.begin(); unique_alleles_it != unique_alleles.end(); ++unique_alleles_it) {
		free(*unique_alleles_it);
	}
}

void Harmonizer2::open_vcf_file(const char* file_name) throw (Harmonizer2Exception) {
	try {
		if (file_name == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_vcf_file( const char* )", __LINE__, 0, "file_name");
		}

		if (strlen(file_name) <= 0) {
			throw Harmonizer2Exception("Harmonizer2", "open_vcf_file( const char* )", __LINE__, 1, "file_name");
		}

		if ((vcf_file != NULL) || (vcf_reader != NULL)) {
			close_vcf_file();
		}

		vcf_file = (char*)malloc((strlen(file_name) + 1u) * sizeof(char));
		if (vcf_file == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_vcf_file( const char* )", __LINE__, 2, (strlen(file_name) + 1u) * sizeof(char));
		}
		strcpy(vcf_file, file_name);

		vcf_reader = ReaderFactory::create(vcf_file);
		vcf_reader->set_file_name(vcf_file);
		vcf_reader->open();
	} catch (ReaderException &e) {
		Harmonizer2Exception new_e(e);
		new_e.add_message("Harmonizer2", "open_vcf_file( const char* )", __LINE__, 3, (vcf_file != NULL) ? vcf_file : "NULL");
		throw new_e;
	} catch (Harmonizer2Exception &e) {
		e.add_message("Harmonizer2", "open_vcf_file( const char* )", __LINE__, 3, (vcf_file != NULL) ? vcf_file : "NULL");
		throw;
	}
}

void Harmonizer2::close_vcf_file() throw (Harmonizer2Exception) {
	try {
		if ((vcf_reader != NULL) && (vcf_reader->is_open())) {
			vcf_reader->close();

			delete vcf_reader;
			vcf_reader = NULL;
		}

		if (vcf_file != NULL) {
			free(vcf_file);
			vcf_file = NULL;
		}

		vcf_file_line_number = 0u;
		vcf_file_column_number = 0;
	} catch (ReaderException &e) {
		Harmonizer2Exception new_e(e);
		new_e.add_message("Harmonizer2", "close_vcf_file()", __LINE__, 4, (vcf_file != NULL) ? vcf_file : "NULL");
		throw new_e;
	}
}


void Harmonizer2::process_vcf_file_header() throw (Harmonizer2Exception) {
	char* line = NULL;
	char* token = NULL;
	int line_length = 0;

	if ((vcf_file == NULL) || (vcf_reader == NULL)) {
		return;
	}

	vcf_file_line_number = 0u;
	vcf_file_column_number = 0;

	try {
		/* Read the first required line with file format description. */
		if ((line_length = vcf_reader->read_line()) > 0) {
			++vcf_file_line_number;
			line = *(vcf_reader->line);
			if ((token = auxiliary::strtok(&line, '=')) != NULL) {
				auxiliary::trim_end(token);
				if (auxiliary::strcmp_ignore_case(token, VCF_FILE_FORMAT) != 0) {
					throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_header()", __LINE__, 8);
				}
			} else {
				throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_header()", __LINE__, 8);
			}
		}

		if (line_length == 0) {
			throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_header()", __LINE__, 9);
		}

		/* Read the mandatory header. Meta-info lines are optional. */
		while ((line_length = vcf_reader->read_line()) > 0) {
			++vcf_file_line_number;
			line = *(vcf_reader->line);

			if (line_length > 1) {
				if (line[0u] != '#') {
					throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_header()", __LINE__, 10);
				}

				if (line[1u] != '#') {
					while ((token = auxiliary::strtok(&line, VCF_FIELD_SEPARATOR)) != NULL) {
						if (vcf_file_column_number < VCF_MANDATORY_COLUMNS_SIZE) {
							if (auxiliary::strcmp_ignore_case(token, vcf_mandatory_columns[vcf_file_column_number]) != 0) {
								throw Harmonizer2Exception("Harmonizer", "process_map_file_header()", __LINE__, 11, vcf_mandatory_columns[vcf_file_column_number], vcf_file_column_number + 1);
							}
						} else {
							/* sample columns */
						}
						++vcf_file_column_number;
					}
					break;
				} else {
					/* process meta-info line if necessary */
				}
			} else {
				throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_header()", __LINE__, 12, vcf_file_line_number);
			}
		}

		if (line_length == 0) {
			throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_header()", __LINE__, 13, vcf_file_line_number + 1u);
		}
	} catch (ReaderException &e) {
		Harmonizer2Exception new_e(e);
		new_e.add_message("Harmonizer2", "process_vcf_file_header()", __LINE__, 5, (vcf_file != NULL) ? vcf_file : "NULL");
		throw new_e;
	}
}

void Harmonizer2::process_vcf_file_data() throw (Harmonizer2Exception) {
	char* line = NULL;
	char* token = NULL;
	int line_length = 0;
	int column_number = 0;

	chr_index* index = NULL;
	position_index_entry* positions_new = NULL;

	char* id = NULL;
	char* chromosome = NULL;
	unsigned long int position = NULL;
	char* ref_allele = NULL;
	char* nonref_allele = NULL;
	unsigned int ref_allele_length = 0u;
	unsigned int nonref_allele_length = 0u;
	char* allele = NULL;
	char variant_type = '\0';
	char variant_subtype = '\0';

	bool vt_found = false;
	int vt_length = strlen(VCF_VARIANT_TYPE);

	if ((vcf_file == NULL) || (vcf_reader == NULL)) {
		return;
	}

	try {
		/* Read data. */
		tokens = (char**)malloc(VCF_MANDATORY_COLUMNS_SIZE * sizeof(char*));
		if (tokens == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_data()", __LINE__, 2, VCF_MANDATORY_COLUMNS_SIZE * sizeof(char*));
		}

		while ((line_length = vcf_reader->read_line()) > 0) {
			++vcf_file_line_number;
			line = *(vcf_reader->line);
			column_number = 0;

			while ((token = auxiliary::strtok(&line, VCF_FIELD_SEPARATOR)) != NULL) {
				if (column_number < VCF_MANDATORY_COLUMNS_SIZE) {
					tokens[column_number] = token;
				}
				++column_number;
			}

			if (column_number < vcf_file_column_number) {
				throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_data()", __LINE__, 14, vcf_file_line_number, vcf_file, column_number, vcf_file_column_number);
			} if (column_number > vcf_file_column_number) {
				throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_data()", __LINE__, 15, vcf_file_line_number, vcf_file, column_number, vcf_file_column_number);
			}

			/* tokens[6] -- filter field. check if all filteres are passed. Must contain "PASS" or "." (if no filters were applied). */
			if ((auxiliary::strcmp_ignore_case(tokens[6u], VCF_PASS) != 0) && (auxiliary::strcmp_ignore_case(tokens[6u], VCF_MISSING) != 0)) {
				continue;
			}

			/* tokens[3] -- ref allele; tokens[4] --check nonref allele. check if polymorphic SNP or small insertion/deletion. */
			ref_allele = tokens[3u];
			ref_allele_length = strlen(ref_allele);
			if ((ref_allele_length <= 0u) || (strspn(ref_allele, "ACGT") != ref_allele_length)) {
				continue;
			}

			nonref_allele = tokens[4u];
			nonref_allele_length = strlen(nonref_allele);
			if (((nonref_allele_length <= 0u) || (strspn(nonref_allele, "ACGT") != nonref_allele_length)) && (auxiliary::strcmp_ignore_case(nonref_allele, VCF_ALT_ALLELE_DEL) != 0)) {
				continue;
			}

			/* tokens[0] -- chromosome. */
			index_by_chr_it = index_by_chr.find(tokens[0u]);
			if (index_by_chr_it != index_by_chr.end()) {
				index = index_by_chr_it->second;
			} else {
				chromosome = (char*)malloc((strlen(tokens[0u]) + 1u) * sizeof(char));
				if (chromosome == NULL) {
					throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_data()", __LINE__, 2, (strlen(tokens[0u]) + 1u) * sizeof(char));
				}
				strcpy(chromosome, tokens[0u]);

				index = new chr_index();

				index->n = 0u;
				index->heap_n = INDEX_HEAP_SIZE;
				index->positions = (position_index_entry*)malloc(index->heap_n * sizeof(position_index_entry));
				if (index->positions == NULL) {
					throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_data()", __LINE__, 2, index->heap_n * sizeof(position_index_entry));
				}

				index_by_chr.insert(pair<char*, chr_index*>(chromosome, index));
			}

			/* tokens[1] -- position. */
			if (!auxiliary::to_ulong_int(tokens[1u], &position)) {
				throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_data()", __LINE__, 20, tokens[1u], vcf_file_line_number, vcf_file);
			}

			/* tokens[2] -- identifier. */
			id = (char*)malloc((strlen(tokens[2u]) + 1u) * sizeof(char));
			if (id == NULL) {
				throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_data()", __LINE__, 2, (strlen(tokens[2u]) + 1u) * sizeof(char));
			}
			strcpy(id, tokens[2u]);

			/* BEGIN: determine variant type (I(indel) and S(snp) are supported) */
			if ((ref_allele_length == 1u) && (nonref_allele_length == 1u)) {
				variant_type = 'S';
				variant_subtype = 'S';
			} else if (auxiliary::strcmp_ignore_case(nonref_allele, VCF_ALT_ALLELE_DEL) == 0) {
				variant_type = 'I';
				variant_subtype = 'D';
			} else if (ref_allele_length > nonref_allele_length) {
				variant_type = 'I';
				variant_subtype = 'D';
			} else if (ref_allele_length < nonref_allele_length) {
				variant_type = 'I';
				variant_subtype = 'I';
			} else {
				continue;
			}
			/* END: determine variant type (I(indel) and S(snp) are supported) */

//			OLD LOGIC
//			/* BEGIN: determine variant type (only I(indel) and S(snp) are supported) */
//			variant_type = 'I';
//
//			/* tokens[7] -- info field that may contain variant type. */
//			while ((token = auxiliary::strtok(&(tokens[7u]), VCF_INFO_FIELD_SEPARATOR)) != NULL) {
//				auxiliary::trim_start(token);
//				if (auxiliary::strcmp_ignore_case(token, VCF_VARIANT_TYPE, vt_length) == 0) {
//					token = strchr(token, '=');
//					if (token != NULL) {
//						++token;
//
//						auxiliary::trim_start(token);
//						auxiliary::trim_end(token);
//
//						if (auxiliary::strcmp_ignore_case(token, VCF_SNP_TYPE) == 0) {
//							variant_type = 'S';
//						}
//					}
//					vt_found = true;
//					break;
//				}
//			}
//
//			/* if info field doesn't contain variant type, then determine it from alleles. */
//			if (!vt_found) {
//				if ((ref_allele_length == 1u) && (nonref_allele_length == 1u)) {
//					variant_type = 'S';
//				}
//			} else {
//				vt_found = false;
//			}
//
//			/* in all other cases it is I (indel) */
//			/* END: determine variant type (only I and S are supported) */

			if (index->n >= index->heap_n) {
				index->heap_n += INDEX_HEAP_INCREMENT;
				positions_new = (position_index_entry*)realloc(index->positions, index->heap_n * sizeof(position_index_entry));
				if (positions_new == NULL) {
					throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_data()", __LINE__, 19, index->heap_n * sizeof(position_index_entry));
				}
				index->positions = positions_new;
				positions_new = NULL;
			}

			index->positions[index->n].position = position;
			index->positions[index->n].id = id;
			index->positions[index->n].type = variant_type;
			index->positions[index->n].subtype = variant_subtype;

			unique_alleles_it = unique_alleles.find(ref_allele);
			if (unique_alleles_it != unique_alleles.end()) {
				allele = *unique_alleles_it;
			} else {
				allele = (char*)malloc((ref_allele_length + 1u) * sizeof(char));
				if (allele == NULL) {
					throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_data()", __LINE__, 2, ((ref_allele_length + 1u) * sizeof(char)));
				}
				strcpy(allele, ref_allele);
				unique_alleles.insert(allele);
			}
			index->positions[index->n].ref_allele = allele;

			unique_alleles_it = unique_alleles.find(nonref_allele);
			if (unique_alleles_it != unique_alleles.end()) {
				allele = *unique_alleles_it;
			} else {
				allele = (char*)malloc((nonref_allele_length + 1u) * sizeof(char));
				if (allele == NULL) {
					throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_data()", __LINE__, 2, ((nonref_allele_length + 1u) * sizeof(char)));
				}
				strcpy(allele, nonref_allele);
				unique_alleles.insert(allele);
			}
			index->positions[index->n].nonref_allele = allele;

			++(index->n);
		}

		if (line_length == 0) {
			throw Harmonizer2Exception("Harmonizer2", "process_vcf_file_data()", __LINE__, 13, vcf_file_line_number + 1);
		}

		free(tokens);
		tokens = NULL;

		index_by_chr_it = index_by_chr.begin();
		while (index_by_chr_it != index_by_chr.end()) {
			index = index_by_chr_it->second;

			qsort(index->positions, index->n, sizeof(position_index_entry), qsort_position_index_entry_cmp);

			++index_by_chr_it;
		}
	} catch (ReaderException &e) {
		Harmonizer2Exception new_e(e);
		new_e.add_message("Harmonizer2", "process_vcf_file_data()", __LINE__, 5, (vcf_file != NULL) ? vcf_file : "NULL");
		throw new_e;
	}
}

void Harmonizer2::index_vcf(const char* file_name) throw (Harmonizer2Exception) {
	try {
		if (file_name == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "index_vcf( const char* )", __LINE__, 0, "file_name");
		}

		if (strlen(file_name) <= 0) {
			throw Harmonizer2Exception("Harmonizer2", "index_vcf( const char* )", __LINE__, 1, "file_name");
		}

		for (index_by_chr_it = index_by_chr.begin(); index_by_chr_it != index_by_chr.end(); ++index_by_chr_it) {
			free(index_by_chr_it->first);
			delete index_by_chr_it->second;
		}

		for (unique_alleles_it = unique_alleles.begin(); unique_alleles_it != unique_alleles.end(); ++unique_alleles_it) {
			free(*unique_alleles_it);
		}

		open_vcf_file(file_name);
		process_vcf_file_header();
		process_vcf_file_data();
		close_vcf_file();
	} catch (Harmonizer2Exception &e) {
		e.add_message("Harmonizer2", "index_vcf( const char* )", __LINE__, 7);
		throw;
	}
}

void Harmonizer2::open_input_file(const char* file_name, const char* id_column_name, const char* chr_column_name, const char* pos_column_name, const char* first_allele_column_name, const char* second_allele_column_name, char field_separator) throw (Harmonizer2Exception) {
	try {
		if ((input_file != NULL) || (input_reader != NULL)) {
			close_input_file();
		}

		if (file_name == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 0, "file_name");
		}

		if (strlen(file_name) <= 0) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 1, "file_name");
		}

		input_file = (char*)malloc((strlen(file_name) + 1u) * sizeof(char));
		if (input_file == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 2, (strlen(file_name) + 1u) * sizeof(char));
		}
		strcpy(input_file, file_name);

		if (id_column_name == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 0, "id_column_name");
		}

		if (strlen(id_column_name) <= 0) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 1, "id_column_name");
		}

		id_column = (char*)malloc((strlen(id_column_name) + 1u) * sizeof(char));
		if (id_column == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 2, (strlen(id_column_name) + 1u) * sizeof(char));
		}
		strcpy(id_column, id_column_name);

		if (chr_column_name == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 0, "chr_column_name");
		}

		if (strlen(chr_column_name) <= 0) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 1, "chr_column_name");
		}

		chr_column = (char*)malloc((strlen(chr_column_name) + 1u) * sizeof(char));
		if (chr_column == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 2, (strlen(chr_column_name) + 1u) * sizeof(char));
		}
		strcpy(chr_column, chr_column_name);

		if (pos_column_name == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 0, "pos_column_name");
		}

		if (strlen(pos_column_name) <= 0) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 1, "pos_column_name");
		}

		pos_column = (char*)malloc((strlen(pos_column_name) + 1u) * sizeof(char));
		if (pos_column == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 2, (strlen(pos_column_name) + 1u) * sizeof(char));
		}
		strcpy(pos_column, pos_column_name);

		if (first_allele_column_name == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 0, "first_allele_column_name");
		}

		if (strlen(first_allele_column_name) <= 0) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 1, "first_allele_column_name");
		}

		first_allele_column = (char*)malloc((strlen(first_allele_column_name) + 1u) * sizeof(char));
		if (first_allele_column == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 2, (strlen(first_allele_column_name) + 1u) * sizeof(char));
		}
		strcpy(first_allele_column, first_allele_column_name);

		if (second_allele_column_name == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 0, "second_allele_column_name");
		}

		if (strlen(second_allele_column_name) <= 0) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 1, "second_allele_column_name");
		}

		second_allele_column = (char*)malloc((strlen(second_allele_column_name) + 1u) * sizeof(char));
		if (second_allele_column == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 2, (strlen(second_allele_column_name) + 1u) * sizeof(char));
		}
		strcpy(second_allele_column, second_allele_column_name);

		separator = field_separator;

		input_reader = ReaderFactory::create(input_file);
		input_reader->set_file_name(input_file);
		input_reader->open();
	} catch (ReaderException &e) {
		Harmonizer2Exception new_e(e);
		new_e.add_message("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 3, (input_file != NULL) ? input_file : "NULL");
		throw new_e;
	} catch (Harmonizer2Exception &e) {
		e.add_message("Harmonizer2", "open_input_file( const char*, const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 3, (input_file != NULL) ? input_file : "NULL");
		throw;
	}
}

void Harmonizer2::open_output_file(const char* file_name, bool gzip) throw (Harmonizer2Exception) {
	try {
		if ((output_file != NULL) || (output_writer != NULL)) {
			close_output_file();
		}

		if (file_name == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_output_file( const char*, bool )", __LINE__, 0, "file_name");
		}

		if (strlen(file_name) <= 0) {
			throw Harmonizer2Exception("Harmonizer2", "open_output_file( const char*, bool )", __LINE__, 1, "file_name");
		}

		output_file = (char*)malloc((strlen(file_name) + 1u) * sizeof(char));
		if (output_file == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_output_file( const char*, bool )", __LINE__, 2, (strlen(file_name) + 1u) * sizeof(char));
		}
		strcpy(output_file, file_name);

		output_writer = WriterFactory::create(gzip ? WriterFactory::GZIP : WriterFactory::TEXT);
		output_writer->set_file_name(output_file);
		output_writer->open();
	} catch (WriterException &e) {
		Harmonizer2Exception new_e(e);
		new_e.add_message("Harmonizer2", "open_output_file( const char*, bool )", __LINE__, 3, (output_file != NULL) ? output_file : "NULL");
		throw new_e;
	} catch (Harmonizer2Exception &e) {
		e.add_message("Harmonizer2", "open_output_file( const char*, bool )", __LINE__, 3, (output_file != NULL) ? output_file : "NULL");
		throw;
	}
}

void Harmonizer2::open_log_file(const char* file_name, bool gzip) throw (Harmonizer2Exception) {
	try {
		if ((log_file != NULL) || (log_writer != NULL)) {
			close_log_file();
		}

		if (file_name == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_log_file( const char*, bool )", __LINE__, 0, "file_name");
		}

		if (strlen(file_name) <= 0) {
			throw Harmonizer2Exception("Harmonizer2", "open_log_file( const char*, bool )", __LINE__, 1, "file_name");
		}

		if (gzip) {
			auxiliary::transform_file_name(&log_file, NULL, file_name, ".log.gz", true);
		} else {
			auxiliary::transform_file_name(&log_file, NULL, file_name, ".log", true);
		}

		if (log_file == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "open_log_file( const char*, bool )", __LINE__, 18);
		}

		log_writer = WriterFactory::create(gzip ? WriterFactory::GZIP : WriterFactory::TEXT);
		log_writer->set_file_name(log_file);
		log_writer->open();
	} catch (WriterException &e) {
		Harmonizer2Exception new_e(e);
		new_e.add_message("Harmonizer2", "open_log_file( const char*, bool )", __LINE__, 3, (log_file != NULL) ? log_file : "NULL");
		throw new_e;
	} catch (Harmonizer2Exception &e) {
		e.add_message("Harmonizer2", "open_log_file( const char*, bool )", __LINE__, 3, (log_file != NULL) ? log_file : "NULL");
		throw;
	}
}

void Harmonizer2::close_input_file() throw (Harmonizer2Exception) {
	try {
		if ((input_reader != NULL) && (input_reader->is_open())) {
			input_reader->close();

			delete input_reader;
			input_reader = NULL;
		}

		if (input_file != NULL) {
			free(input_file);
			input_file = NULL;
		}

		if (chr_column != NULL) {
			free(chr_column);
			chr_column = NULL;
		}

		if (id_column != NULL) {
			free(id_column);
			id_column = NULL;
		}

		if (first_allele_column != NULL) {
			free(first_allele_column);
			first_allele_column = NULL;
		}

		if (second_allele_column != NULL) {
			free(second_allele_column);
			second_allele_column = NULL;
		}

		separator = '\0';

		if (header_backup != NULL) {
			free(header_backup);
			header_backup = NULL;
		}

		input_file_line_number = 0u;
		input_file_column_number = 0;

		id_column_pos = numeric_limits<int>::min();
		second_allele_column_pos = numeric_limits<int>::min();
		first_allele_column_pos = numeric_limits<int>::min();
	} catch (ReaderException &e) {
		Harmonizer2Exception new_e(e);
		new_e.add_message("Harmonizer2", "close_input_file()", __LINE__, 4, (input_file != NULL) ? input_file : "NULL");
		throw new_e;
	}
}

void Harmonizer2::close_output_file() throw (Harmonizer2Exception) {
	try {
		if (output_writer != NULL) {
			output_writer->close();

			delete output_writer;
			output_writer = NULL;
		}

		if (output_file != NULL) {
			free(output_file);
			output_file = NULL;
		}
	} catch (WriterException &e) {
		Harmonizer2Exception new_e(e);
		new_e.add_message("Harmonizer2", "close_output_file()", __LINE__, 4, (output_file != NULL) ? output_file : "NULL");
		throw new_e;
	}
}

void Harmonizer2::close_log_file() throw (Harmonizer2Exception) {
	try {
		if (log_writer != NULL) {
			log_writer->close();

			delete log_writer;
			log_writer = NULL;
		}

		if (log_file != NULL) {
			free(log_file);
			log_file = NULL;
		}
	} catch (WriterException &e) {
		Harmonizer2Exception new_e(e);
		new_e.add_message("Harmonizer2", "close_log_file()", __LINE__, 4, (log_file != NULL) ? log_file : "NULL");
		throw new_e;
	}
}

void Harmonizer2::process_header() throw (Harmonizer2Exception) {
	char* header = NULL;
	char* token = NULL;

	input_file_line_number = 0u;
	input_file_column_number = 0;

	id_column_pos = numeric_limits<int>::min();
	chr_column_pos = numeric_limits<int>::min();
	pos_column_pos = numeric_limits<int>::min();
	first_allele_column_pos = numeric_limits<int>::min();
	second_allele_column_pos = numeric_limits<int>::min();

	if ((input_file == NULL) || (input_reader == NULL) ||
			(id_column == NULL) || (chr_column == NULL) || (pos_column == NULL) ||
			(first_allele_column == NULL) || (second_allele_column == NULL)) {
		return;
	}

	try {
		if (input_reader->read_line() <= 0) {
			throw Harmonizer2Exception("Harmonizer2", "process_header()", __LINE__, 16, 1, input_file);
		}

		header = *(input_reader->line);

		header_backup = (char*)malloc((strlen(header) + 1u) * sizeof(char));
		if (header_backup == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "process_header()", __LINE__, 2, ((strlen(header) + 1u) * sizeof(char)));
		}
		strcpy(header_backup, header);

		while ((token = auxiliary::strtok(&header, separator)) != NULL) {
			if (auxiliary::strcmp_ignore_case(token, id_column) == 0) {
				id_column_pos = input_file_column_number;
			} else if (auxiliary::strcmp_ignore_case(token, chr_column) == 0) {
				chr_column_pos = input_file_column_number;
			} else if (auxiliary::strcmp_ignore_case(token, pos_column) == 0) {
				pos_column_pos = input_file_column_number;
			} else if (auxiliary::strcmp_ignore_case(token, first_allele_column) == 0) {
				first_allele_column_pos = input_file_column_number;
			} else if (auxiliary::strcmp_ignore_case(token, second_allele_column) == 0) {
				second_allele_column_pos = input_file_column_number;
			}
			++input_file_column_number;
		}

		if (id_column_pos < 0) {
			throw Harmonizer2Exception("Harmonizer2", "process_header()", __LINE__, 17, id_column, input_file);
		}

		if (chr_column_pos < 0) {
			throw Harmonizer2Exception("Harmonizer2", "process_header()", __LINE__, 17, chr_column, input_file);
		}

		if (pos_column_pos < 0) {
			throw Harmonizer2Exception("Harmonizer2", "process_header()", __LINE__, 17, pos_column, input_file);
		}

		if (first_allele_column_pos < 0) {
			throw Harmonizer2Exception("Harmonizer2", "process_header()", __LINE__, 17, first_allele_column, input_file);
		}

		if (second_allele_column_pos < 0) {
			throw Harmonizer2Exception("Harmonizer2", "process_header()", __LINE__, 17, second_allele_column, input_file);
		}

		++input_file_line_number;
	} catch (ReaderException &e) {
		Harmonizer2Exception new_e(e);
		new_e.add_message("Harmonizer2", "process_header()", __LINE__, 5, (input_file != NULL) ? input_file : "NULL");
		throw new_e;
	}
}

void Harmonizer2::write_columns() throw (WriterException) {
	output_writer->write("%s", tokens[0]);
	for (int column_number = 1; column_number < input_file_column_number; ++column_number) {
		output_writer->write("%c%s", separator, tokens[column_number]);
	}
	output_writer->write("\n");
}

void Harmonizer2::harmonize(bool flip, bool drop) throw (Harmonizer2Exception) {
	char* line = NULL;
	char* token = NULL;
	int line_length = 0;
	int column_number = 0;

	char* buffer = NULL;

	chr_index* index = NULL;
	unsigned int index_size = 0u;

	char* end_ptr = NULL;

	char* first_allele = NULL;
	char* second_allele = NULL;
	unsigned int first_allele_length = 0u;
	unsigned int second_allele_length = 0u;
	char c_first_allele = '\0';
	char c_second_allele = '\0';
	bool empty_first_allele = false;
	bool empty_second_allele = false;
	bool recode_alleles = true;
	char flipped_first_allele[2] = {'\0', '\0'};
	char flipped_second_allele[2] = {'\0', '\0'};

	position_index_entry* positions = NULL;
	position_index_entry query_position;
	position_index_entry* found_position = NULL;

	long int start = 0;
	long int end = 0;

	bool type_ok = false;
	bool alleles_ok = false;
	bool strand_ok = true;

	unsigned int ref_allele_length = 0u;
	unsigned int nonref_allele_length = 0u;

	if ((input_file == NULL) || (input_reader == NULL) ||
			(output_file == NULL) || (output_writer == NULL) ||
			(log_file == NULL) || (log_writer == NULL) ||
			(id_column == NULL) || (id_column_pos < 0) ||
			(chr_column == NULL) || (chr_column_pos < 0) ||
			(pos_column == NULL) || (pos_column_pos < 0) ||
			(first_allele_column == NULL) || (first_allele_column_pos < 0) ||
			(second_allele_column == NULL) || (second_allele_column_pos < 0)) {
		return;
	}

	try {
		tokens = (char**)malloc(input_file_column_number * sizeof(char*));
		if (tokens == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "harmonize( bool )", __LINE__, 2, input_file_column_number * sizeof(char*));
		}

		buffer = (char*)malloc(16777216 * sizeof(char));
		if (buffer == NULL) {
			throw Harmonizer2Exception("Harmonizer2", "harmonize( bool )", __LINE__, 2, 16777216 * sizeof(char));
		}

		output_writer->write("%s\n", header_backup);

		while ((line_length = input_reader->read_line()) > 0) {
			line = *(input_reader->line);
			++input_file_line_number;

			column_number = 0;
			while ((token = auxiliary::strtok(&line, separator)) != NULL) {
				if (column_number < input_file_column_number) {
					tokens[column_number] = token;
				}
				++column_number;
			}

			if (column_number < input_file_column_number) {
				throw Harmonizer2Exception("Harmonizer2", "harmonize( bool )", __LINE__, 14, input_file_line_number, input_file, column_number, input_file_column_number);
			} else if (column_number > input_file_column_number) {
				throw Harmonizer2Exception("Harmonizer2", "harmonize( bool )", __LINE__, 15, input_file_line_number, input_file, column_number, input_file_column_number);
			}

			auxiliary::trim(&(tokens[chr_column_pos]));
			index_by_chr_it = index_by_chr.find(tokens[chr_column_pos]);
			if (index_by_chr_it == index_by_chr.end()) {
				log_writer->write("Line %u: (WARNING) Chromosome of %s at %s:%s is not in VCF.\n", input_file_line_number, tokens[id_column_pos], tokens[chr_column_pos], tokens[pos_column_pos]);
				if (!drop) write_columns();
				continue;
			}

			index = index_by_chr_it->second;
			positions = index->positions;
			index_size = index->n;

			query_position.position = strtoul(tokens[pos_column_pos], &end_ptr, 10);
			if (end_ptr == tokens[pos_column_pos]) {
				log_writer->write("Line %u: (WARNING) Position %s of %s can't be parsed to integer.\n", input_file_line_number, tokens[pos_column_pos], tokens[id_column_pos]);
				if (!drop) write_columns();
				continue;
			}

			query_position.subtype = '\0';

			auxiliary::trim(&(tokens[first_allele_column_pos]));
			first_allele = tokens[first_allele_column_pos];
			first_allele_length = strlen(first_allele);

			auxiliary::trim(&(tokens[second_allele_column_pos]));
			second_allele = tokens[second_allele_column_pos];
			second_allele_length = strlen(second_allele);

			c_first_allele = '\0';
			c_second_allele = '\0';
			empty_first_allele = false;
			empty_second_allele = false;
			recode_alleles = true;

			if ((first_allele_length == 1) && (second_allele_length == 1)) {
				c_first_allele = toupper(first_allele[0]);
				c_second_allele = toupper(second_allele[0]);

				if (((c_first_allele == 'A') || (c_first_allele == 'C') ||	(c_first_allele == 'G') || (c_first_allele == 'T')) &&
						((c_second_allele == 'A') || (c_second_allele == 'C') || (c_second_allele == 'G') || (c_second_allele == 'T'))) {
					query_position.type = 'S';
					query_position.subtype = 'S';
				} else if (((c_first_allele == 'D') && (c_second_allele == 'R')) ||	((c_first_allele == 'R') && (c_second_allele == 'D'))) {
					query_position.type = 'I';
					query_position.subtype = 'D';
					recode_alleles = false;
				} else if (((c_first_allele == 'I') && (c_second_allele == 'R')) ||	((c_first_allele == 'R') && (c_second_allele ==  'I'))) {
					query_position.type = 'I';
					query_position.subtype = 'I';
					recode_alleles = false;
				} else if (((c_first_allele == '.') || (c_first_allele == '-')) &&
						((c_second_allele == 'A') || (c_second_allele == 'C') || (c_second_allele == 'G') || (c_second_allele == 'T'))) {
					query_position.type = 'I';
					empty_first_allele = true;
				} else if (((c_second_allele == '.') || (c_second_allele == '-')) &&
						((c_first_allele == 'A') || (c_first_allele == 'C') || (c_first_allele == 'G') || (c_first_allele == 'T'))) {
					query_position.type = 'I';
					empty_second_allele = true;
				} else {
					log_writer->write("Line %u: (WARNING) Alleles %s/%s of %s at %s:%s don't match any variation type.\n", input_file_line_number, tokens[first_allele_column_pos], tokens[second_allele_column_pos], tokens[id_column_pos], tokens[chr_column_pos], tokens[pos_column_pos]);
					if (!drop) write_columns();
					continue;
				}
			} else {
				query_position.type = 'I';

				if ((auxiliary::strcmp_ignore_case(first_allele, ".") == 0) || (auxiliary::strcmp_ignore_case(first_allele, "-") == 0) ||
						(auxiliary::strcmp_ignore_case(first_allele, "NA") == 0) || (auxiliary::strcmp_ignore_case(first_allele, VCF_ALT_ALLELE_DEL) == 0)) {
					empty_first_allele = true;
				} else if (strspn(first_allele, "ACGT") != first_allele_length) {
					log_writer->write("Line %u: (WARNING) Alleles %s/%s of %s at %s:%s don't match any variation type.\n", input_file_line_number, tokens[first_allele_column_pos], tokens[second_allele_column_pos], tokens[id_column_pos], tokens[chr_column_pos], tokens[pos_column_pos]);
					if (!drop) write_columns();
					continue;
				}

				if ((auxiliary::strcmp_ignore_case(second_allele, ".") == 0) || (auxiliary::strcmp_ignore_case(second_allele, "-") == 0) ||
						(auxiliary::strcmp_ignore_case(second_allele, "NA") == 0) || (auxiliary::strcmp_ignore_case(second_allele, VCF_ALT_ALLELE_DEL) == 0)) {
					empty_second_allele = true;
				} else if (strspn(second_allele, "ACGT") != second_allele_length) {
					log_writer->write("Line %u: (WARNING) Alleles %s/%s of %s at %s:%s don't match any variation type.\n", input_file_line_number, tokens[first_allele_column_pos], tokens[second_allele_column_pos], tokens[id_column_pos], tokens[chr_column_pos], tokens[pos_column_pos]);
					if (!drop) write_columns();
					continue;
				}

				if (empty_first_allele && empty_second_allele) {
					log_writer->write("Line %u: (WARNING) Alleles %s/%s of %s at %s:%s don't match any variation type.\n", input_file_line_number, tokens[first_allele_column_pos], tokens[second_allele_column_pos], tokens[id_column_pos], tokens[chr_column_pos], tokens[pos_column_pos]);
					if (!drop) write_columns();
					continue;
				}

				if (!empty_first_allele && !empty_second_allele && (first_allele_length == second_allele_length)) {
					log_writer->write("Line %u: (WARNING) Alleles %s/%s of %s at %s:%s don't match any variation type.\n", input_file_line_number, tokens[first_allele_column_pos], tokens[second_allele_column_pos], tokens[id_column_pos], tokens[chr_column_pos], tokens[pos_column_pos]);
					if (!drop) write_columns();
					continue;
				}
			}

			found_position = (position_index_entry*)bsearch(&query_position, positions, index_size, sizeof(position_index_entry), qsort_position_index_entry_cmp);
			if (found_position == NULL) {
				log_writer->write("Line %u: (WARNING) Position %s:%s of %s is not in VCF.\n", input_file_line_number, tokens[chr_column_pos], tokens[pos_column_pos], tokens[id_column_pos]);
				if (!drop) write_columns();
				continue;
			}

			start = found_position - positions;
			end = start;

			while ((start >= 0) && (qsort_position_index_entry_cmp(&(positions[start]), &query_position) == 0)) {
				--start;
			}
			++start;

			while ((++end < index_size) && (qsort_position_index_entry_cmp(&(positions[end]), &query_position) == 0));
			--end;

			type_ok = false;
			alleles_ok = false;
			strand_ok = true;

			while (start <= end) {
				if (positions[start].type == query_position.type) {
					type_ok = true;

					if (query_position.type == 'S') { // SNPs
						if (((auxiliary::strcmp_ignore_case(positions[start].ref_allele, first_allele) == 0) && (auxiliary::strcmp_ignore_case(positions[start].nonref_allele, second_allele) == 0)) ||
								((auxiliary::strcmp_ignore_case(positions[start].ref_allele, second_allele) == 0) && (auxiliary::strcmp_ignore_case(positions[start].nonref_allele, first_allele) == 0))) {
							alleles_ok = true;
							break;
						}

						// Check for strand mismatch
						if (c_first_allele == 'A') {
							flipped_first_allele[0] = 'T';
						} else if (c_first_allele == 'C') {
							flipped_first_allele[0] = 'G';
						} else if (c_first_allele == 'G') {
							flipped_first_allele[0] = 'C';
						} else {
							flipped_first_allele[0] = 'A';
						}

						if (c_second_allele == 'A') {
							flipped_second_allele[0] = 'T';
						} else if (c_second_allele == 'C') {
							flipped_second_allele[0] = 'G';
						} else if (c_second_allele == 'G') {
							flipped_second_allele[0] = 'C';
						} else {
							flipped_second_allele[0] = 'A';
						}

						if (((auxiliary::strcmp_ignore_case(positions[start].ref_allele, flipped_first_allele) == 0) && (auxiliary::strcmp_ignore_case(positions[start].nonref_allele, flipped_second_allele) == 0)) ||
								((auxiliary::strcmp_ignore_case(positions[start].ref_allele, flipped_second_allele) == 0) && (auxiliary::strcmp_ignore_case(positions[start].nonref_allele, flipped_first_allele) == 0))) {
							if (flip) {
								log_writer->write("Line %u: Strand of %s at %s:%s flipped.\n", input_file_line_number, tokens[id_column_pos], tokens[chr_column_pos], tokens[pos_column_pos]);

								tokens[first_allele_column_pos][0] = flipped_first_allele[0];
								tokens[first_allele_column_pos][1] = '\0';
								tokens[second_allele_column_pos][0] = flipped_second_allele[0];
								tokens[second_allele_column_pos][1] = '\0';
							} else {
								strand_ok = false;
							}
							alleles_ok = true;
							break;
						}
					} else { // INDELs
						if (recode_alleles) { // Check for alleles match and recode to R/I/D
							if (auxiliary::strcmp_ignore_case(positions[start].nonref_allele, VCF_ALT_ALLELE_DEL) == 0) { // Structural variant (i.e. <DEL>) in alternate VCF allele.
								if (empty_second_allele && (auxiliary::strcmp_ignore_case(positions[start].ref_allele, first_allele) == 0)) {
									log_writer->write("Line %u: Alleles %s/%s of %s at %s:%s changed to %s/%s.\n", input_file_line_number, tokens[first_allele_column_pos], tokens[second_allele_column_pos], tokens[id_column_pos], tokens[chr_column_pos], tokens[pos_column_pos], "R", "D");

									tokens[first_allele_column_pos][0] = 'R';
									tokens[first_allele_column_pos][1] = '\0';
									tokens[second_allele_column_pos][0] = 'D';
									tokens[second_allele_column_pos][1] = '\0';

									alleles_ok = true;
									break;
								} else if (empty_first_allele && (auxiliary::strcmp_ignore_case(positions[start].ref_allele, second_allele) == 0)) {
									log_writer->write("Line %u: Alleles %s/%s of %s at %s:%s changed to %s/%s.\n", input_file_line_number, tokens[first_allele_column_pos], tokens[second_allele_column_pos], tokens[id_column_pos], tokens[chr_column_pos], tokens[pos_column_pos], "D", "R");

									tokens[first_allele_column_pos][0] = 'D';
									tokens[first_allele_column_pos][1] = '\0';
									tokens[second_allele_column_pos][0] = 'R';
									tokens[second_allele_column_pos][1] = '\0';

									alleles_ok = true;
									break;
								}
							} else {
								if ((auxiliary::strcmp_ignore_case(positions[start].ref_allele, first_allele) == 0) && (auxiliary::strcmp_ignore_case(positions[start].nonref_allele, second_allele) == 0)) {
									if (positions[start].subtype == 'D') {
										log_writer->write("Line %u: Alleles %s/%s of %s at %s:%s changed to %s/%s.\n", input_file_line_number, tokens[first_allele_column_pos], tokens[second_allele_column_pos], tokens[id_column_pos], tokens[chr_column_pos], tokens[pos_column_pos], "R", "D");

										tokens[first_allele_column_pos][0] = 'R';
										tokens[first_allele_column_pos][1] = '\0';
										tokens[second_allele_column_pos][0] = 'D';
										tokens[second_allele_column_pos][1] = '\0';
									} else {
										log_writer->write("Line %u: Alleles %s/%s of %s at %s:%s changed to %s/%s.\n", input_file_line_number, tokens[first_allele_column_pos], tokens[second_allele_column_pos], tokens[id_column_pos], tokens[chr_column_pos], tokens[pos_column_pos], "R", "I");

										tokens[first_allele_column_pos][0] = 'R';
										tokens[first_allele_column_pos][1] = '\0';
										tokens[second_allele_column_pos][0] = 'I';
										tokens[second_allele_column_pos][1] = '\0';
									}

									alleles_ok = true;
									break;
								} else if ((auxiliary::strcmp_ignore_case(positions[start].ref_allele, second_allele) == 0) && (auxiliary::strcmp_ignore_case(positions[start].nonref_allele, first_allele) == 0)) {
									if (positions[start].subtype == 'D') {
										log_writer->write("Line %u: Alleles %s/%s of %s at %s:%s changed to %s/%s.\n", input_file_line_number, tokens[first_allele_column_pos], tokens[second_allele_column_pos], tokens[id_column_pos], tokens[chr_column_pos], tokens[pos_column_pos], "D", "R");

										tokens[first_allele_column_pos][0] = 'D';
										tokens[first_allele_column_pos][1] = '\0';
										tokens[second_allele_column_pos][0] = 'R';
										tokens[second_allele_column_pos][1] = '\0';
									} else {
										log_writer->write("Line %u: Alleles %s/%s of %s at %s:%s changed to %s/%s.\n", input_file_line_number, tokens[first_allele_column_pos], tokens[second_allele_column_pos], tokens[id_column_pos], tokens[chr_column_pos], tokens[pos_column_pos], "I", "R");

										tokens[first_allele_column_pos][0] = 'I';
										tokens[first_allele_column_pos][1] = '\0';
										tokens[second_allele_column_pos][0] = 'R';
										tokens[second_allele_column_pos][1] = '\0';
									}

									alleles_ok = true;
									break;
								}
							}
						} else { // Already coded as R/I/D. Only check for sub-type (insertion/deletion) consistency.
							if (positions[start].subtype == query_position.subtype) {
								alleles_ok = true;
								break;
							}
						}
					}
				}

				++start;
			}

			if (!type_ok) {
				log_writer->write("Line %u: (WARNING) Type of %s at %s:%s doesn't match type in VCF.\n", input_file_line_number, tokens[id_column_pos], tokens[chr_column_pos], tokens[pos_column_pos]);
				if (!drop) write_columns();
				continue;
			} else if (!alleles_ok) {
				log_writer->write("Line %u: (WARNING) Alleles %s/%s of %s at %s:%s don't match alleles in VCF.\n", input_file_line_number, tokens[first_allele_column_pos], tokens[second_allele_column_pos], tokens[id_column_pos], tokens[chr_column_pos], tokens[pos_column_pos]);
				if (!drop) write_columns();
				continue;
			} else if (!strand_ok) {
				log_writer->write("Line %u: (WARNING) Strand of %s at %s:%s doesn't match strand in VCF.\n", input_file_line_number, tokens[id_column_pos], tokens[chr_column_pos], tokens[pos_column_pos]);
				if (!drop) write_columns();
				continue;
			}

			found_position = &(positions[start]);
			if (found_position->type == 'S') {
				sprintf(buffer, "%s:%lu:%s", tokens[chr_column_pos], found_position->position, VCF_SNP_TYPE);
			} else {
				sprintf(buffer, "%s:%lu:%s", tokens[chr_column_pos], found_position->position, VCF_INDEL_TYPE_01);
			}

			if (strcmp(buffer, tokens[id_column_pos]) != 0) {
				log_writer->write("Line %u: %s changed to %s.\n", input_file_line_number, tokens[id_column_pos], buffer);
				tokens[id_column_pos] = buffer;
			}

			write_columns();
		}

		if (line_length == 0) {
			throw Harmonizer2Exception("Harmonizer2", "harmonize( bool )", __LINE__, 21, input_file_line_number + 1, input_file);
		}

		free(tokens);
		tokens = NULL;

		free(buffer);
		buffer = NULL;
	} catch (ReaderException &e) {
		Harmonizer2Exception new_e(e);
		new_e.add_message("Harmonizer2", "harmonize( bool )", __LINE__, 5, (input_file != NULL) ? input_file : "NULL");
		throw new_e;
	} catch (WriterException &e) {
		Harmonizer2Exception new_e(e);
		new_e.add_message("Harmonizer2", "harmonize( bool )", __LINE__, 5, (input_file != NULL) ? input_file : "NULL");
		throw new_e;
	} catch (Harmonizer2Exception &e) {
		e.add_message("Harmonizer2", "harmonize( bool )", __LINE__, 5, (input_file != NULL) ? input_file : "NULL");
		throw;
	}
}

inline int Harmonizer2::qsort_position_index_entry_cmp(const void* first, const void* second) {
	position_index_entry* first_entry = (position_index_entry*)first;
	position_index_entry* second_entry = (position_index_entry*)second;

	if (first_entry->position > second_entry->position) {
		return 1;
	} else if (first_entry->position < second_entry->position) {
		return -1;
	}

	return 0;
}
