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

#include "include/Harmonizer.h"

const char Harmonizer::VCF_FIELD_SEPARATOR = '\t';
const char* Harmonizer::VCF_FILE_FORMAT = "##fileformat";
const char* Harmonizer::VCF_CHROM = "#CHROM";
const char* Harmonizer::VCF_POS = "POS";
const char* Harmonizer::VCF_ID = "ID";
const char* Harmonizer::VCF_REF = "REF";
const char* Harmonizer::VCF_ALT = "ALT";
const char* Harmonizer::VCF_QUAL = "QUAL";
const char* Harmonizer::VCF_FILTER = "FILTER";
const char* Harmonizer::VCF_INFO = "INFO";
const char* Harmonizer::VCF_FORMAT = "FORMAT";
const int Harmonizer::VCF_MANDATORY_COLUMNS_SIZE = 9;
const char* Harmonizer::vcf_mandatory_columns[VCF_MANDATORY_COLUMNS_SIZE] = {
		VCF_CHROM, VCF_POS, VCF_ID, VCF_REF, VCF_ALT, VCF_QUAL, VCF_FILTER, VCF_INFO, VCF_FORMAT
};
const char* Harmonizer::VCF_PASS = "PASS";
const char* Harmonizer::VCF_MISSING = ".";
const char Harmonizer::VCF_INFO_FIELD_SEPARATOR = ';';
const char* Harmonizer::VCF_VARIANT_TYPE = "VT";
const char* Harmonizer::VCF_SNP_TYPE = "SNP";
const char* Harmonizer::VCF_INDEL_TYPE_01 = "INDEL";
const char* Harmonizer::VCF_INDEL_TYPE_02 = "I";
const char* Harmonizer::VCF_ALT_ALLELE_DEL = "<DEL>";
const char* Harmonizer::VCF_ALT_ALLELE_INS = "<INS>";

const unsigned int Harmonizer::MAP_HEAP_SIZE = 3000000;
const unsigned int Harmonizer::MAP_HEAP_INCREMENT = 1000000;

Harmonizer::Harmonizer() : map_file(NULL), map_reader(NULL), map_file_line_number(0u), map_file_column_number(0),
		input_file(NULL), output_file(NULL), log_file(NULL), reader(NULL), writer(NULL), log_writer(NULL), tokens(NULL),
		chr_column(NULL), id_column(NULL), ref_allele_column(NULL), nonref_allele_column(NULL),
		separator('\0'), header_backup(NULL), file_column_number(0),
		chr_column_pos(numeric_limits<int>::min()), id_column_pos(numeric_limits<int>::min()),
		ref_allele_column_pos(numeric_limits<int>::min()), nonref_allele_column_pos(numeric_limits<int>::min()),
		map_index_by_chr(auxiliary::bool_strcmp_ignore_case) {

}

Harmonizer::~Harmonizer() {
	if (map_file != NULL) {
		free(map_file);
		map_file = NULL;
	}

	if (map_reader != NULL) {
		delete map_reader;
		map_reader = NULL;
	}

	if (input_file != NULL) {
		free(input_file);
		input_file = NULL;
	}

	if (output_file != NULL) {
		free(output_file);
		output_file = NULL;
	}

	if (log_file != NULL) {
		free(log_file);
		log_file = NULL;
	}

	if (reader != NULL) {
		delete reader;
		reader = NULL;
	}

	if (writer != NULL) {
		delete writer;
		writer = NULL;
	}

	if (log_writer != NULL) {
		delete log_writer;
		log_writer = NULL;
	}

	if (tokens != NULL) {
		delete tokens;
		tokens = NULL;
	}

	if (chr_column != NULL) {
		free(chr_column);
		chr_column = NULL;
	}

	if (id_column != NULL) {
		free(id_column);
		id_column = NULL;
	}

	if (ref_allele_column != NULL) {
		free(ref_allele_column);
		ref_allele_column = NULL;
	}

	if (nonref_allele_column != NULL) {
		free(nonref_allele_column);
		nonref_allele_column = NULL;
	}

	if (header_backup != NULL) {
		free(header_backup);
		header_backup = NULL;
	}

	for (map_index_by_chr_it = map_index_by_chr.begin(); map_index_by_chr_it != map_index_by_chr.end(); ++map_index_by_chr_it) {
		free(map_index_by_chr_it->first);
		delete map_index_by_chr_it->second;
	}
}

void Harmonizer::open_map_file(const char* file_name) throw (HarmonizerException) {
	try {
		if (file_name == NULL) {
			throw HarmonizerException("Harmonizer", "open_map_file( const char* )", __LINE__, 0, "file_name");
		}

		if (strlen(file_name) <= 0) {
			throw HarmonizerException("Harmonizer", "open_map_file( const char* )", __LINE__, 1, "file_name");
		}

		if ((map_file != NULL) || (map_reader != NULL)) {
			close_map_file();
		}

		map_file = (char*)malloc((strlen(file_name) + 1u) * sizeof(char));
		if (map_file == NULL) {
			throw HarmonizerException("Harmonizer", "open_map_file( const char* )", __LINE__, 2, (strlen(file_name) + 1u) * sizeof(char));
		}
		strcpy(map_file, file_name);

		map_reader = ReaderFactory::create(map_file);
		map_reader->set_file_name(map_file);
		map_reader->open();
	} catch (ReaderException &e) {
		HarmonizerException new_e(e);
		new_e.add_message("Harmonizer", "open_map_file( const char* )", __LINE__, 3, (map_file != NULL) ? map_file : "NULL");
		throw new_e;
	} catch (HarmonizerException &e) {
		e.add_message("Harmonizer", "open_map_file( const char* )", __LINE__, 3, (map_file != NULL) ? map_file : "NULL");
		throw;
	}
}

void Harmonizer::close_map_file() throw (HarmonizerException) {
	try {
		if ((map_reader != NULL) && (map_reader->is_open())) {
			map_reader->close();

			delete map_reader;
			map_reader = NULL;
		}

		if (map_file != NULL) {
			free(map_file);
			map_file = NULL;
		}

		map_file_line_number = 0u;
		map_file_column_number = 0;
	} catch (ReaderException &e) {
		HarmonizerException new_e(e);
		new_e.add_message("Harmonizer", "close_map_file()", __LINE__, 4, (map_file != NULL) ? map_file : "NULL");
		throw new_e;
	}
}


void Harmonizer::process_map_file_header() throw (HarmonizerException) {
	char* line = NULL;
	int line_length = 0;

	char* token = NULL;

	if ((map_file == NULL) || (map_reader == NULL)) {
		return;
	}

	map_file_line_number = 0u;
	map_file_column_number = 0;

	try {
		/* Read the first required line with file format description. */
		if ((line_length = map_reader->read_line()) > 0) {
			++map_file_line_number;
			line = *(map_reader->line);
			if ((token = auxiliary::strtok(&line, '=')) != NULL) {
				auxiliary::trim_end(token);
				if (auxiliary::strcmp_ignore_case(token, VCF_FILE_FORMAT) != 0) {
					throw HarmonizerException("Harmonizer", "process_map_file_header()", __LINE__, 8);
				}
			} else {
				throw HarmonizerException("Harmonizer", "process_map_file_header()", __LINE__, 8);
			}
		}

		if (line_length == 0) {
			throw HarmonizerException("Harmonizer", "process_map_file_header()", __LINE__, 9);
		}

		/* Read the mandatory header. Meta-info lines are optional. */
		while ((line_length = map_reader->read_line()) > 0) {
			++map_file_line_number;
			line = *(map_reader->line);

			if (line_length > 1) {
				if (line[0u] != '#') {
					throw HarmonizerException("Harmonizer", "process_map_file_header()", __LINE__, 10);
				}

				if (line[1u] != '#') {
					while ((token = auxiliary::strtok(&line, VCF_FIELD_SEPARATOR)) != NULL) {
						if (map_file_column_number < VCF_MANDATORY_COLUMNS_SIZE) {
							if (auxiliary::strcmp_ignore_case(token, vcf_mandatory_columns[map_file_column_number]) != 0) {
								throw HarmonizerException("Harmonizer", "process_map_file_header()", __LINE__, 11, vcf_mandatory_columns[map_file_column_number], map_file_column_number + 1u);
							}
						} else {
							/* sample columns */
						}
						++map_file_column_number;
					}
					break;
				} else {
					/* process meta-info line if necessary */

//					if (auxiliary::strcmp_ignore_case(line, "##ALT", 5) == 0) {
//						token = strstr(line, "ID");
//						if (token == NULL) {
//							continue;
//						}
//
//						if (strpbrk(token, ",>") != NULL) {
//							strpbrk(token, ",>")[0u] = '\0';
//						}
//
//						token = strchr(token, '=');
//						if (token == NULL) {
//							continue;
//						}
//
//						++token;
//
//						auxiliary::trim(&token);
//						if (strlen(token) <= 0) {
//							continue;
//						}
//
//						/* the token contains variant type name */
//					}
				}
			} else {
				throw HarmonizerException("Harmonizer", "process_map_file_header()", __LINE__, 12, map_file_line_number);
			}
		}

		if (line_length == 0) {
			throw HarmonizerException("Harmonizer", "process_map_file_header()", __LINE__, 13, map_file_line_number + 1u);
		}
	} catch (ReaderException &e) {
		HarmonizerException new_e(e);
		new_e.add_message("Harmonizer", "process_map_file_header()", __LINE__, 5, (map_file != NULL) ? map_file : "NULL");
		throw new_e;
	}
}

void Harmonizer::process_map_file_data() throw (HarmonizerException) {
	char* line = NULL;
	int line_length = 0;

	int column_number = 0;

	char* token = NULL;

	if ((map_file == NULL) || (map_reader == NULL)) {
		return;
	}

	chr_index* index = NULL;
	position_index_entry* positions_new = NULL;

	char* chromosome = NULL;
	char* id = NULL;
	unsigned long int position = NULL;
	char type = '\0';
	char* ref_allele = NULL;
	unsigned int ref_allele_length = 0u;
	char* nonref_allele = NULL;
	unsigned int nonref_allele_length = 0u;

	bool vt_found = false;
	int vt_length = strlen(VCF_VARIANT_TYPE);

	try {
		/* Read data. */
		tokens = (char**)malloc(VCF_MANDATORY_COLUMNS_SIZE * sizeof(char*));
		if (tokens == NULL) {
			throw HarmonizerException("Harmonizer", "process_map_file_data()", __LINE__, 2, VCF_MANDATORY_COLUMNS_SIZE * sizeof(char*));
		}

		while ((line_length = map_reader->read_line()) > 0) {
			++map_file_line_number;
			line = *(map_reader->line);
			column_number = 0;

			while ((token = auxiliary::strtok(&line, VCF_FIELD_SEPARATOR)) != NULL) {
				if (column_number < VCF_MANDATORY_COLUMNS_SIZE) {
					tokens[column_number] = token;
				}
				++column_number;
			}

			if (column_number < map_file_column_number) {
				throw HarmonizerException("Harmonizer", "process_map_file_data()", __LINE__, 14, map_file_line_number, map_file, column_number, map_file_column_number);
			} if (column_number > map_file_column_number) {
				throw HarmonizerException("Harmonizer", "process_map_file_data()", __LINE__, 15, map_file_line_number, map_file, column_number, map_file_column_number);
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
			if (((nonref_allele_length <= 0u) || (strspn(nonref_allele, "ACGT") != nonref_allele_length)) &&
					(auxiliary::strcmp_ignore_case(nonref_allele, VCF_ALT_ALLELE_DEL) != 0) &&
					(auxiliary::strcmp_ignore_case(nonref_allele, VCF_ALT_ALLELE_INS) != 0)) {
				continue;
			}

			/* tokens[0] -- chromosome. */
			map_index_by_chr_it = map_index_by_chr.find(tokens[0u]);
			if (map_index_by_chr_it != map_index_by_chr.end()) {
				index = map_index_by_chr_it->second;
			} else {
				chromosome = (char*)malloc((strlen(tokens[0u]) + 1u) * sizeof(char));
				if (chromosome == NULL) {
					throw HarmonizerException("Harmonizer", "process_map_file_data()", __LINE__, 2, (strlen(tokens[0u]) + 1u) * sizeof(char));
				}
				strcpy(chromosome, tokens[0u]);

				index = new chr_index();

				index->n = 0u;
				index->heap_n = MAP_HEAP_SIZE;
				index->positions = (position_index_entry*)malloc(index->heap_n * sizeof(position_index_entry));
				if (index->positions == NULL) {
					throw HarmonizerException("Harmonizer", "process_map_file_data()", __LINE__, 2, index->heap_n * sizeof(position_index_entry));
				}
				index->ids = NULL;

				map_index_by_chr.insert(pair<char*, chr_index*>(chromosome, index));
			}

			/* tokens[1] -- position. */
			if (!auxiliary::to_ulong_int(tokens[1u], &position)) {
				throw HarmonizerException("Harmonizer", "process_map_file_data()", __LINE__, 20, tokens[1u], map_file_line_number, map_file);
			}

			/* tokens[2] -- identifier. */
			id = (char*)malloc((strlen(tokens[2u]) + 1u) * sizeof(char));
			if (id == NULL) {
				throw HarmonizerException("Harmonizer", "process_map_file_data()", __LINE__, 2, (strlen(tokens[2u]) + 1u) * sizeof(char));
			}
			strcpy(id, tokens[2u]);

			/* BEGIN: determine variant type (only I(indel) and S(snp) are supported) */
			type = 'I';

			/* tokens[7] -- info field that may contain variant type. */
			while ((token = auxiliary::strtok(&(tokens[7u]), VCF_INFO_FIELD_SEPARATOR)) != NULL) {
				auxiliary::trim_start(token);
				if (auxiliary::strcmp_ignore_case(token, VCF_VARIANT_TYPE, vt_length) == 0) {
					token = strchr(token, '=');
					if (token != NULL) {
						++token;

						auxiliary::trim_start(token);
						auxiliary::trim_end(token);

						if (auxiliary::strcmp_ignore_case(token, VCF_SNP_TYPE) == 0) {
							type = 'S';
						}
					}
					vt_found = true;
					break;
				}
			}

			/* if info field doesn't contain variant type, then determine it from alleles. */
			if (!vt_found) {
				if ((ref_allele_length == 1u) && (nonref_allele_length == 1u)) {
					type = 'S';
				}
			} else {
				vt_found = false;
			}

			/* in all other cases it is I (indel) */
			/* END: determine variant type (only I and S are supported) */

			if (index->n >= index->heap_n) {
				index->heap_n += MAP_HEAP_INCREMENT;
				positions_new = (position_index_entry*)realloc(index->positions, index->heap_n * sizeof(position_index_entry));
				if (positions_new == NULL) {
					throw HarmonizerException("Harmonizer", "process_map_file_data()", __LINE__, 19, index->heap_n * sizeof(position_index_entry));
				}
				index->positions = positions_new;
				positions_new = NULL;
			}

			index->positions[index->n].position = position;
			index->positions[index->n].id = id;
			index->positions[index->n].type = type;
			if (type == 'I') {
				index->positions[index->n].ref_allele = (char*)malloc((ref_allele_length + 1u) * sizeof(char));
				if (index->positions[index->n].ref_allele == NULL) {
					throw HarmonizerException("Harmonizer", "process_map_file_data()", __LINE__, 2, ((ref_allele_length + 1u) * sizeof(char)));
				}
				strcpy(index->positions[index->n].ref_allele, ref_allele);

				index->positions[index->n].nonref_allele = (char*)malloc((nonref_allele_length + 1u) * sizeof(char));
				if (index->positions[index->n].nonref_allele == NULL) {
					throw HarmonizerException("Harmonizer", "process_map_file_data()", __LINE__, 2, ((nonref_allele_length + 1u) * sizeof(char)));
				}
				strcpy(index->positions[index->n].nonref_allele, nonref_allele);
			} else {
				index->positions[index->n].ref_allele = NULL;
				index->positions[index->n].nonref_allele = NULL;
			}

			++(index->n);
		}

		if (line_length == 0) {
			throw HarmonizerException("Harmonizer", "process_map_file_data()", __LINE__, 13, map_file_line_number + 1u);
		}

		free(tokens);
		tokens = NULL;

		map_index_by_chr_it = map_index_by_chr.begin();
		while (map_index_by_chr_it != map_index_by_chr.end()) {
			index = map_index_by_chr_it->second;

			qsort(index->positions, index->n, sizeof(position_index_entry), qsort_position_index_entry_cmp);

			index->ids = (id_index_entry*)malloc(index->n * sizeof(id_index_entry));
			if (index->ids == NULL) {
				throw HarmonizerException("Harmonizer", "process_map_file_data()", __LINE__, 2, index->n * sizeof(id_index_entry));
			}

			for (unsigned int i = 0u; i < index->n; ++i) {
				index->ids[i].id = index->positions[i].id;
				index->ids[i].location = i;
			}

			qsort(index->ids, index->n, sizeof(id_index_entry), qsort_id_index_entry_cmp);

			++map_index_by_chr_it;
		}
	} catch (ReaderException &e) {
		HarmonizerException new_e(e);
		new_e.add_message("Harmonizer", "process_map_file_data()", __LINE__, 5, (map_file != NULL) ? map_file : "NULL");
		throw new_e;
	}
}

void Harmonizer::index_map(const char* file_name) throw (HarmonizerException) {
	try {
		if (file_name == NULL) {
			throw HarmonizerException("Harmonizer", "index_map( const char* )", __LINE__, 0, "file_name");
		}

		if (strlen(file_name) <= 0) {
			throw HarmonizerException("Harmonizer", "index_map( const char* )", __LINE__, 1, "file_name");
		}

		for (map_index_by_chr_it = map_index_by_chr.begin(); map_index_by_chr_it != map_index_by_chr.end(); ++map_index_by_chr_it) {
			free(map_index_by_chr_it->first);
			delete map_index_by_chr_it->second;
		}

		open_map_file(file_name);
		process_map_file_header();
		process_map_file_data();
		close_map_file();
	} catch (HarmonizerException &e) {
		e.add_message("Harmonizer", "index_map( const char* )", __LINE__, 7);
		throw;
	}
}

void Harmonizer::open_input_file(const char* file_name, const char* chr_column_name, const char* id_column_name, const char* ref_allele_column_name, const char* nonref_allele_column_name, char field_separator) throw (HarmonizerException) {
	try {
		if ((input_file != NULL) || (reader != NULL)) {
			close_input_file();
		}

		if (file_name == NULL) {
			throw HarmonizerException("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 0, "file_name");
		}

		if (strlen(file_name) <= 0) {
			throw HarmonizerException("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 1, "file_name");
		}

		input_file = (char*)malloc((strlen(file_name) + 1u) * sizeof(char));
		if (input_file == NULL) {
			throw HarmonizerException("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 2, (strlen(file_name) + 1u) * sizeof(char));
		}
		strcpy(input_file, file_name);

		if (chr_column_name == NULL) {
			throw HarmonizerException("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 0, "chr_column_name");
		}

		if (strlen(chr_column_name) <= 0) {
			throw HarmonizerException("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 1, "chr_column_name");
		}

		chr_column = (char*)malloc((strlen(chr_column_name) + 1u) * sizeof(char));
		if (chr_column == NULL) {
			throw HarmonizerException("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 2, (strlen(chr_column_name) + 1u) * sizeof(char));
		}
		strcpy(chr_column, chr_column_name);

		if (id_column_name == NULL) {
			throw HarmonizerException("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 0, "id_column_name");
		}

		if (strlen(id_column_name) <= 0) {
			throw HarmonizerException("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 1, "id_column_name");
		}

		id_column = (char*)malloc((strlen(id_column_name) + 1u) * sizeof(char));
		if (id_column == NULL) {
			throw HarmonizerException("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 2, (strlen(id_column_name) + 1u) * sizeof(char));
		}
		strcpy(id_column, id_column_name);

		if (ref_allele_column_name == NULL) {
			throw HarmonizerException("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 0, "ref_allele_column_name");
		}

		if (strlen(ref_allele_column_name) <= 0) {
			throw HarmonizerException("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 1, "ref_allele_column_name");
		}

		ref_allele_column = (char*)malloc((strlen(ref_allele_column_name) + 1u) * sizeof(char));
		if (ref_allele_column == NULL) {
			throw HarmonizerException("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 2, (strlen(ref_allele_column_name) + 1u) * sizeof(char));
		}
		strcpy(ref_allele_column, ref_allele_column_name);

		if (nonref_allele_column_name == NULL) {
			throw HarmonizerException("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 0, "nonref_allele_column_name");
		}

		if (strlen(nonref_allele_column_name) <= 0) {
			throw HarmonizerException("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 1, "nonref_allele_column_name");
		}

		nonref_allele_column = (char*)malloc((strlen(nonref_allele_column_name) + 1u) * sizeof(char));
		if (nonref_allele_column == NULL) {
			throw HarmonizerException("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 2, (strlen(nonref_allele_column_name) + 1u) * sizeof(char));
		}
		strcpy(nonref_allele_column, nonref_allele_column_name);

		separator = field_separator;

		reader = ReaderFactory::create(input_file);
		reader->set_file_name(input_file);
		reader->open();
	} catch (ReaderException &e) {
		HarmonizerException new_e(e);
		new_e.add_message("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 3, (input_file != NULL) ? input_file : "NULL");
		throw new_e;
	} catch (HarmonizerException &e) {
		e.add_message("Harmonizer", "open_input_file( const char*, const char*, const char*, const char*, const char*, char )", __LINE__, 3, (input_file != NULL) ? input_file : "NULL");
		throw;
	}
}

void Harmonizer::open_output_file(const char* file_name, bool gzip) throw (HarmonizerException) {
	try {
		if ((output_file != NULL) || (writer != NULL)) {
			close_output_file();
		}

		if (file_name == NULL) {
			throw HarmonizerException("Harmonizer", "open_output_file( const char*, bool )", __LINE__, 0, "file_name");
		}

		if (strlen(file_name) <= 0) {
			throw HarmonizerException("Harmonizer", "open_output_file( const char*, bool )", __LINE__, 1, "file_name");
		}

		output_file = (char*)malloc((strlen(file_name) + 1u) * sizeof(char));
		if (output_file == NULL) {
			throw HarmonizerException("Harmonizer", "open_output_file( const char*, bool )", __LINE__, 2, (strlen(file_name) + 1u) * sizeof(char));
		}
		strcpy(output_file, file_name);

		writer = WriterFactory::create(gzip ? WriterFactory::GZIP : WriterFactory::TEXT);
		writer->set_file_name(output_file);
		writer->open();
	} catch (WriterException &e) {
		HarmonizerException new_e(e);
		new_e.add_message("Harmonizer", "open_output_file( const char*, bool )", __LINE__, 3, (output_file != NULL) ? output_file : "NULL");
		throw new_e;
	} catch (HarmonizerException &e) {
		e.add_message("Harmonizer", "open_output_file( const char*, bool )", __LINE__, 3, (output_file != NULL) ? output_file : "NULL");
		throw;
	}
}

void Harmonizer::open_log_file(const char* file_name, bool gzip) throw (HarmonizerException) {
	try {
		if ((log_file != NULL) || (log_writer != NULL)) {
			close_log_file();
		}

		if (file_name == NULL) {
			throw HarmonizerException("Harmonizer", "open_log_file( const char*, bool )", __LINE__, 0, "file_name");
		}

		if (strlen(file_name) <= 0) {
			throw HarmonizerException("Harmonizer", "open_log_file( const char*, bool )", __LINE__, 1, "file_name");
		}

		if (gzip) {
			auxiliary::transform_file_name(&log_file, NULL, file_name, ".log.gz", true);
		} else {
			auxiliary::transform_file_name(&log_file, NULL, file_name, ".log", true);
		}

		if (log_file == NULL) {
			throw HarmonizerException("Harmonizer", "open_log_file( const char*, bool )", __LINE__, 18);
		}

		log_writer = WriterFactory::create(gzip ? WriterFactory::GZIP : WriterFactory::TEXT);
		log_writer->set_file_name(log_file);
		log_writer->open();
	} catch (WriterException &e) {
		HarmonizerException new_e(e);
		new_e.add_message("Harmonizer", "open_log_file( const char*, bool )", __LINE__, 3, (log_file != NULL) ? log_file : "NULL");
		throw new_e;
	} catch (HarmonizerException &e) {
		e.add_message("Harmonizer", "open_log_file( const char*, bool )", __LINE__, 3, (log_file != NULL) ? log_file : "NULL");
		throw;
	}
}

void Harmonizer::close_input_file() throw (HarmonizerException) {
	try {
		if ((reader != NULL) && (reader->is_open())) {
			reader->close();

			delete reader;
			reader = NULL;
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

		if (ref_allele_column != NULL) {
			free(ref_allele_column);
			ref_allele_column = NULL;
		}

		if (nonref_allele_column != NULL) {
			free(nonref_allele_column);
			nonref_allele_column = NULL;
		}

		separator = '\0';

		if (header_backup != NULL) {
			free(header_backup);
			header_backup = NULL;
		}

		file_column_number = 0;
		id_column_pos = numeric_limits<int>::min();
		ref_allele_column_pos = numeric_limits<int>::min();
		nonref_allele_column_pos = numeric_limits<int>::min();
	} catch (ReaderException &e) {
		HarmonizerException new_e(e);
		new_e.add_message("Harmonizer", "close_input_file()", __LINE__, 4, (input_file != NULL) ? input_file : "NULL");
		throw new_e;
	}
}

void Harmonizer::close_output_file() throw (HarmonizerException) {
	try {
		if (writer != NULL) {
			writer->close();

			delete writer;
			writer = NULL;
		}

		if (output_file != NULL) {
			free(output_file);
			output_file = NULL;
		}
	} catch (WriterException &e) {
		HarmonizerException new_e(e);
		new_e.add_message("Harmonizer", "close_output_file()", __LINE__, 4, (output_file != NULL) ? output_file : "NULL");
		throw new_e;
	}
}

void Harmonizer::close_log_file() throw (HarmonizerException) {
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
		HarmonizerException new_e(e);
		new_e.add_message("Harmonizer", "close_log_file()", __LINE__, 4, (log_file != NULL) ? log_file : "NULL");
		throw new_e;
	}
}

void Harmonizer::process_header() throw (HarmonizerException) {
	char* header = NULL;
	char* token = NULL;

	if ((input_file == NULL) || (reader == NULL) || (chr_column == NULL) || (id_column == NULL)) {
		return;
	}

	file_column_number = 0;
	chr_column_pos = numeric_limits<int>::min();
	id_column_pos = numeric_limits<int>::min();
	ref_allele_column_pos = numeric_limits<int>::min();
	nonref_allele_column_pos = numeric_limits<int>::min();

	try {
		if (reader->read_line() <= 0) {
			throw HarmonizerException("Harmonizer", "process_header()", __LINE__, 16, 1, input_file);
		}

		header = *(reader->line);

		header_backup = (char*)malloc((strlen(header) + 1u) * sizeof(char));
		if (header_backup == NULL) {
			throw HarmonizerException("Harmonizer", "process_header()", __LINE__, 2, ((strlen(header) + 1u) * sizeof(char)));
		}
		strcpy(header_backup, header);

		while ((token = auxiliary::strtok(&header, separator)) != NULL) {
			if (auxiliary::strcmp_ignore_case(token, chr_column) == 0) {
				chr_column_pos = file_column_number;
			} else if (auxiliary::strcmp_ignore_case(token, id_column) == 0) {
				id_column_pos = file_column_number;
			} else if (auxiliary::strcmp_ignore_case(token, ref_allele_column) == 0) {
				ref_allele_column_pos = file_column_number;
			} else if (auxiliary::strcmp_ignore_case(token, nonref_allele_column) == 0) {
				nonref_allele_column_pos = file_column_number;
			}
			++file_column_number;
		}

		if (ref_allele_column_pos < 0) {
			throw HarmonizerException("Harmonizer", "process_header()", __LINE__, 17, ref_allele_column, input_file);
		}

		if (nonref_allele_column_pos < 0) {
			throw HarmonizerException("Harmonizer", "process_header()", __LINE__, 17, nonref_allele_column, input_file);
		}


		if (chr_column_pos < 0) {
			throw HarmonizerException("Harmonizer", "process_header()", __LINE__, 17, chr_column, input_file);
		}

		if (id_column_pos < 0) {
			throw HarmonizerException("Harmonizer", "process_header()", __LINE__, 17, id_column, input_file);
		}
	} catch (ReaderException &e) {
		HarmonizerException new_e(e);
		new_e.add_message("Harmonizer", "process_header()", __LINE__, 5, (input_file != NULL) ? input_file : "NULL");
		throw new_e;
	}
}

void Harmonizer::write_columns() throw (WriterException) {
	writer->write("%s", tokens[0]);
	for (int column_number = 1; column_number < file_column_number; ++column_number) {
		writer->write("%c%s", separator, tokens[column_number]);
	}
	writer->write("\n");
}

void Harmonizer::harmonize(bool drop) throw (HarmonizerException) {
	char* line = NULL;
	int line_length = 0;
	unsigned int line_number = 1u;

	char* token = NULL;

	int column_number = 0;

	chr_index* index = NULL;
	unsigned int index_size = 0u;

	position_index_entry* position_index = NULL;
	position_index_entry query_position;
	position_index_entry* found_position = NULL;

	id_index_entry* id_index = NULL;
	id_index_entry query_id;
	id_index_entry* found_id = NULL;

	unsigned int found_id_index_entry_pos = 0u;
	unsigned int found_position_index_entry_pos = 0u;
	long int position = 0;

	char* ref_allele = NULL;
	char* nonref_allele = NULL;
	unsigned int ref_allele_length = 0u;
	unsigned int nonref_allele_length = 0u;

	char* buffer = NULL;
	char* buffer_tmp = NULL;
	unsigned int buffer_tokens_number = 0u;
	char* buffer_token = NULL;
	char* buffer_tokens[3u];
	char* end_ptr = NULL;

	bool empty_ref_allele = false;
	bool empty_nonref_allele = false;
	bool recode_alleles = false;
	bool consistent_type = false;
	bool consistent_alleles = false;

	if ((input_file == NULL) || (reader == NULL) ||
			(output_file == NULL) || (writer == NULL) ||
			(log_file == NULL) || (log_writer == NULL) ||
			(chr_column == NULL) || (chr_column_pos < 0) ||
			(id_column == NULL) || (id_column_pos < 0) ||
			(ref_allele_column == NULL) || (ref_allele_column_pos < 0) ||
			(nonref_allele_column == NULL) || (nonref_allele_column_pos < 0)) {
		return;
	}

	try {
		tokens = (char**)malloc(file_column_number * sizeof(char*));
		if (tokens == NULL) {
			throw HarmonizerException("Harmonizer", "harmonize( bool )", __LINE__, 2, file_column_number * sizeof(char*));
		}

		buffer = (char*)malloc(16777216 * sizeof(char));
		if (buffer == NULL) {
			throw HarmonizerException("Harmonizer", "harmonize( bool )", __LINE__, 2, 16777216 * sizeof(char));
		}

		writer->write("%s\n", header_backup);

		while ((line_length = reader->read_line()) > 0) {
			line = *(reader->line);
			++line_number;

			recode_alleles = false;
			consistent_type = false;
			consistent_alleles = false;

			column_number = 0;
			while ((token = auxiliary::strtok(&line, separator)) != NULL) {
				if (column_number < file_column_number) {
					tokens[column_number] = token;
				}
				++column_number;
			}

			if (column_number < file_column_number) {
				throw HarmonizerException("Harmonizer", "harmonize( bool )", __LINE__, 14, line_number, input_file, column_number, file_column_number);
			} else if (column_number > file_column_number) {
				throw HarmonizerException("Harmonizer", "harmonize( bool )", __LINE__, 15, line_number, input_file, column_number, file_column_number);
			}

			auxiliary::trim(&(tokens[chr_column_pos]));
			map_index_by_chr_it = map_index_by_chr.find(tokens[chr_column_pos]);
			if (map_index_by_chr_it == map_index_by_chr.end()) {
				log_writer->write("Line %u: (WARNING) Chromosome %s of %s is not in the map.\n", line_number, tokens[chr_column_pos], tokens[id_column_pos]);
				if (!drop) write_columns();
				continue;
			}

			index = map_index_by_chr_it->second;
			index_size = index->n;

			position_index = index->positions;
			id_index = index->ids;

			found_id = NULL;
			found_position = NULL;

			auxiliary::trim(&(tokens[ref_allele_column_pos]));
			ref_allele = tokens[ref_allele_column_pos];
			ref_allele_length = strlen(ref_allele);

			auxiliary::trim(&(tokens[nonref_allele_column_pos]));
			nonref_allele = tokens[nonref_allele_column_pos];
			nonref_allele_length = strlen(nonref_allele);

			auxiliary::trim(&(tokens[id_column_pos]));
			if (auxiliary::strcmp_ignore_case(tokens[id_column_pos], "rs", 2) == 0) { // Look up by rsID
				query_id.id = tokens[id_column_pos];
				found_id = (id_index_entry*)bsearch(&query_id, id_index, index_size, sizeof(id_index_entry), qsort_id_index_entry_cmp);
				if (found_id == NULL) {
					log_writer->write("Line %u: (WARNING) %s is not in the map.\n", line_number, query_id.id);
					if (!drop) write_columns();
					continue;
				}

				if (((auxiliary::strcmp_ignore_case(ref_allele, "A") == 0) || (auxiliary::strcmp_ignore_case(ref_allele, "C") == 0) ||
					(auxiliary::strcmp_ignore_case(ref_allele, "G") == 0) || (auxiliary::strcmp_ignore_case(ref_allele, "T") == 0)) &&
					((auxiliary::strcmp_ignore_case(nonref_allele, "A") == 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "C") == 0) ||
						(auxiliary::strcmp_ignore_case(nonref_allele, "G") == 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "T") == 0))) {
					query_id.type = 'S';
				} else if (((auxiliary::strcmp_ignore_case(ref_allele, "D") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "R") == 0)) ||
							((auxiliary::strcmp_ignore_case(ref_allele, "R") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "D") == 0))) {
					query_id.type = 'I';
				} else if (((auxiliary::strcmp_ignore_case(ref_allele, "I") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "R") == 0)) ||
							((auxiliary::strcmp_ignore_case(ref_allele, "R") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "I") == 0))) {
					query_id.type = 'I';
				} else {
					query_id.type = 'I';

					if (auxiliary::strcmp_ignore_case(ref_allele, nonref_allele) == 0) {
						log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}

					if ((ref_allele_length == nonref_allele_length) || (strspn(ref_allele, "ACGT") != ref_allele_length) || (strspn(nonref_allele, "ACGT") != nonref_allele_length)) {
						log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}

					recode_alleles = true;
				}

				found_id_index_entry_pos = found_id - id_index;

				if (recode_alleles) {
					position = found_id_index_entry_pos;
					while ((position >= 0) && (qsort_id_index_entry_cmp(&(id_index[position]), &query_id) == 0)) {
						if (position_index[id_index[position].location].type == query_id.type) {
							consistent_type = true;

							if ((auxiliary::strcmp_ignore_case(position_index[id_index[position].location].ref_allele, ref_allele) == 0) &&
									(auxiliary::strcmp_ignore_case(position_index[id_index[position].location].nonref_allele, nonref_allele) == 0)) {
								consistent_alleles = true;

								if (ref_allele_length > nonref_allele_length) {
									log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "R", "D");
									tokens[ref_allele_column_pos][0] = 'R';
									tokens[ref_allele_column_pos][1] = '\0';
									tokens[nonref_allele_column_pos][0] = 'D';
									tokens[nonref_allele_column_pos][1] = '\0';
								} else if (ref_allele_length < nonref_allele_length) {
									log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "R", "I");
									tokens[ref_allele_column_pos][0] = 'R';
									tokens[ref_allele_column_pos][1] = '\0';
									tokens[nonref_allele_column_pos][0] = 'I';
									tokens[nonref_allele_column_pos][1] = '\0';
								}

								break;
							} else if ((auxiliary::strcmp_ignore_case(position_index[id_index[position].location].ref_allele, nonref_allele) == 0) &&
									(auxiliary::strcmp_ignore_case(position_index[id_index[position].location].nonref_allele, ref_allele) == 0)) {
								consistent_alleles = true;

								if (ref_allele_length > nonref_allele_length) {
									log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "I", "R");
									tokens[ref_allele_column_pos][0] = 'I';
									tokens[ref_allele_column_pos][1] = '\0';
									tokens[nonref_allele_column_pos][0] = 'R';
									tokens[nonref_allele_column_pos][1] = '\0';
								} else if (ref_allele_length < nonref_allele_length) {
									log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "D", "R");
									tokens[ref_allele_column_pos][0] = 'D';
									tokens[ref_allele_column_pos][1] = '\0';
									tokens[nonref_allele_column_pos][0] = 'R';
									tokens[nonref_allele_column_pos][1] = '\0';
								}

								break;
							}
						}
						--position;
					}

					if (!consistent_type) {
						position = found_id_index_entry_pos;
						++position;
						while ((position < index_size) && (qsort_id_index_entry_cmp(&(id_index[position]), &query_id) == 0)) {
							if (position_index[id_index[position].location].type == query_id.type) {
								consistent_type = true;

								if ((auxiliary::strcmp_ignore_case(position_index[id_index[position].location].ref_allele, ref_allele) == 0) &&
										(auxiliary::strcmp_ignore_case(position_index[id_index[position].location].nonref_allele, nonref_allele) == 0)) {
									consistent_alleles = true;

									if (ref_allele_length > nonref_allele_length) {
										log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "R", "D");
										tokens[ref_allele_column_pos][0] = 'R';
										tokens[ref_allele_column_pos][1] = '\0';
										tokens[nonref_allele_column_pos][0] = 'D';
										tokens[nonref_allele_column_pos][1] = '\0';
									} else if (ref_allele_length < nonref_allele_length) {
										log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "R", "I");
										tokens[ref_allele_column_pos][0] = 'R';
										tokens[ref_allele_column_pos][1] = '\0';
										tokens[nonref_allele_column_pos][0] = 'I';
										tokens[nonref_allele_column_pos][1] = '\0';
									}

									break;
								} else if ((auxiliary::strcmp_ignore_case(position_index[id_index[position].location].ref_allele, nonref_allele) == 0) &&
										(auxiliary::strcmp_ignore_case(position_index[id_index[position].location].nonref_allele, ref_allele) == 0)) {
									consistent_alleles = true;

									if (ref_allele_length > nonref_allele_length) {
										log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "I", "R");
										tokens[ref_allele_column_pos][0] = 'I';
										tokens[ref_allele_column_pos][1] = '\0';
										tokens[nonref_allele_column_pos][0] = 'R';
										tokens[nonref_allele_column_pos][1] = '\0';
									} else if (ref_allele_length < nonref_allele_length) {
										log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "D", "R");
										tokens[ref_allele_column_pos][0] = 'D';
										tokens[ref_allele_column_pos][1] = '\0';
										tokens[nonref_allele_column_pos][0] = 'R';
										tokens[nonref_allele_column_pos][1] = '\0';
									}

									break;
								}
							}
							++position;
						}
					}

					if (!consistent_type) {
						log_writer->write("Line %u: (WARNING) Type of %s doesn't match type in the map.\n", line_number, query_id.id);
						if (!drop) write_columns();
						continue;
					}

					if (!consistent_alleles) {
						log_writer->write("Line %u: (WARNING) Alleles of %s don't match alleles in the map.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}
				} else {
					position = found_id_index_entry_pos;
					while ((position >= 0) && (qsort_id_index_entry_cmp(&(id_index[position]), &query_id) == 0)) {
						if (position_index[id_index[position].location].type == query_id.type) {
							consistent_type = true;
							break;
						}
						--position;
					}

					if (!consistent_type) {
						position = found_id_index_entry_pos;
						++position;
						while ((position < index_size) && (qsort_id_index_entry_cmp(&(id_index[position]), &query_id) == 0)) {
							if (position_index[id_index[position].location].type == query_id.type) {
								consistent_type = true;
								break;
							}
							++position;
						}
					}

					if (!consistent_type) {
						log_writer->write("Line %u: (WARNING) Type of %s doesn't match type in the map.\n", line_number, query_id.id);
						if (!drop) write_columns();
						continue;
					}
				}

				found_position = &(position_index[id_index[position].location]);
				if (found_position->type == 'S') {
					sprintf(buffer, "%s:%lu:%s", tokens[chr_column_pos], found_position->position, VCF_SNP_TYPE);
				} else {
					sprintf(buffer, "%s:%lu:%s", tokens[chr_column_pos], found_position->position, VCF_INDEL_TYPE_01);
				}
				log_writer->write("Line %u: %s changed to %s.\n", line_number, tokens[id_column_pos], buffer);
				tokens[id_column_pos] = buffer;
			} else if (auxiliary::strcmp_ignore_case(tokens[id_column_pos], "merged_del", 10) == 0) { // Look up by id
				query_id.id = tokens[id_column_pos];
				found_id = (id_index_entry*)bsearch(&query_id, id_index, index_size, sizeof(id_index_entry), qsort_id_index_entry_cmp);
				if (found_id == NULL) {
					log_writer->write("Line %u: (WARNING) %s is not in the map.\n", line_number, query_id.id);
					if (!drop) write_columns();
					continue;
				}

				if (((auxiliary::strcmp_ignore_case(ref_allele, "D") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "R") == 0)) ||
							((auxiliary::strcmp_ignore_case(ref_allele, "R") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "D") == 0))) {
					query_id.type = 'I';
				} else {
					query_id.type = 'I';

					empty_ref_allele = false;
					empty_nonref_allele = false;

					if ((strcmp(ref_allele, "-") == 0) || (strcmp(ref_allele, ".") == 0) || (auxiliary::strcmp_ignore_case(ref_allele, "NA") == 0)) {
						empty_ref_allele = true;
					} else if (strspn(ref_allele, "ACGT") != ref_allele_length) {
						log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}

					if ((strcmp(nonref_allele, "-") == 0) || (strcmp(nonref_allele, ".") == 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "NA") == 0)) {
						empty_nonref_allele = true;
					} else if (strspn(nonref_allele, "ACGT") != nonref_allele_length) {
						log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}

					if (empty_ref_allele && empty_nonref_allele) {
						log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}

					if (!empty_ref_allele && !empty_nonref_allele && (ref_allele_length == nonref_allele_length)) {
						log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}

					if (empty_ref_allele) {
						log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "D", "R");
						tokens[ref_allele_column_pos][0] = 'D';
						tokens[ref_allele_column_pos][1] = '\0';
						tokens[nonref_allele_column_pos][0] = 'R';
						tokens[nonref_allele_column_pos][1] = '\0';
					} else if (empty_nonref_allele) {
						log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "R", "D");
						tokens[ref_allele_column_pos][0] = 'R';
						tokens[ref_allele_column_pos][1] = '\0';
						tokens[nonref_allele_column_pos][0] = 'D';
						tokens[nonref_allele_column_pos][1] = '\0';
					} else if (ref_allele_length > nonref_allele_length) {
						log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "R", "D");
						tokens[ref_allele_column_pos][0] = 'R';
						tokens[ref_allele_column_pos][1] = '\0';
						tokens[nonref_allele_column_pos][0] = 'D';
						tokens[nonref_allele_column_pos][1] = '\0';
					} else if (ref_allele_length < nonref_allele_length) {
						log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "D", "R");
						tokens[ref_allele_column_pos][0] = 'D';
						tokens[ref_allele_column_pos][1] = '\0';
						tokens[nonref_allele_column_pos][0] = 'R';
						tokens[nonref_allele_column_pos][1] = '\0';
					} else {
						log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}
				}

				found_id_index_entry_pos = found_id - id_index;

				position = found_id_index_entry_pos;
				while ((position >= 0) && (qsort_id_index_entry_cmp(&(id_index[position]), &query_id) == 0)) {
					if (position_index[id_index[position].location].type == query_id.type) {
						consistent_type = true;
						break;
					}
					--position;
				}

				if (!consistent_type) {
					position = found_id_index_entry_pos;
					++position;
					while ((position < index_size) && (qsort_id_index_entry_cmp(&(id_index[position]), &query_id) == 0)) {
						if (position_index[id_index[position].location].type == query_id.type) {
							consistent_type = true;
							break;
						}
						++position;
					}
				}

				if (!consistent_type) {
					log_writer->write("Line %u: (WARNING) Type of %s doesn't match type in the map.\n", line_number, query_id.id);
					if (!drop) write_columns();
					continue;
				}

				found_position = &(position_index[id_index[position].location]);
				if (found_position->type == 'S') {
					sprintf(buffer, "%s:%lu:%s", tokens[chr_column_pos], found_position->position, VCF_SNP_TYPE);
				} else {
					sprintf(buffer, "%s:%lu:%s", tokens[chr_column_pos], found_position->position, VCF_INDEL_TYPE_01);
				}
				log_writer->write("Line %u: %s changed to %s.\n", line_number, tokens[id_column_pos], buffer);
				tokens[id_column_pos] = buffer;
			} else { // Look up by position
				strcpy(buffer, tokens[id_column_pos]);
				buffer_tmp = buffer;

				buffer_tokens_number = 0u;
				buffer_tokens[0u] = buffer_tokens[1u] = buffer_tokens[2u] = NULL;
				while ((buffer_tokens_number < 3u) && ((buffer_token = auxiliary::strtok(&buffer_tmp, ':')) != NULL)) {
					buffer_tokens[buffer_tokens_number++] = buffer_token;
				}

				if ((buffer_tokens_number < 2u) || (buffer_tokens_number > 3u)) {
					log_writer->write("Line %u: (WARNING) %s identifier doesn't follow CHR:POS, CHR:POS:TYPE, or rsID formats.\n", line_number, tokens[id_column_pos]);
					if (!drop) write_columns();
					continue;
				}

				query_position.position = strtoul(buffer_tokens[1u], &end_ptr, 10);
				if (end_ptr == buffer_tokens[1u]) {
					log_writer->write("Line %u: (WARNING) Position of %s can't be parsed to integer.\n", line_number, tokens[id_column_pos]);
					if (!drop) write_columns();
					continue;
				}

				found_position = (position_index_entry*)bsearch(&query_position, position_index, index_size, sizeof(position_index_entry), qsort_position_index_entry_cmp);
				if (found_position == NULL) {
					log_writer->write("Line %u: (WARNING) %s is not in the map.\n", line_number, tokens[id_column_pos]);
					if (!drop) write_columns();
					continue;
				}

				query_position.type = '\0';
				if (buffer_tokens[2u] != NULL) {
					if (auxiliary::strcmp_ignore_case(buffer_tokens[2u], VCF_SNP_TYPE) == 0) {
						//CHR:POS:SNP
						query_position.type = 'S';
						if (((auxiliary::strcmp_ignore_case(ref_allele, "A") != 0) && (auxiliary::strcmp_ignore_case(ref_allele, "C") != 0) &&
								(auxiliary::strcmp_ignore_case(ref_allele, "G") != 0) && (auxiliary::strcmp_ignore_case(ref_allele, "T") != 0)) ||
								((auxiliary::strcmp_ignore_case(nonref_allele, "A") != 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "C") != 0) &&
										(auxiliary::strcmp_ignore_case(nonref_allele, "G") != 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "T") != 0))) {
							log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
							if (!drop) write_columns();
							continue;
						}
					} else if (auxiliary::strcmp_ignore_case(buffer_tokens[2u], "D") == 0) {
						//CHR:POS:D
						query_position.type = 'I';
						if ((strspn(ref_allele, "ACGT") == ref_allele_length) && ((strspn(nonref_allele, "ACGT") == nonref_allele_length))) {
							if (ref_allele_length > nonref_allele_length) {
								log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "R", "D");
								tokens[ref_allele_column_pos][0] = 'R';
								tokens[ref_allele_column_pos][1] = '\0';
								tokens[nonref_allele_column_pos][0] = 'D';
								tokens[nonref_allele_column_pos][1] = '\0';
							} else if (ref_allele_length < nonref_allele_length) {
								log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "D", "R");
								tokens[ref_allele_column_pos][0] = 'D';
								tokens[ref_allele_column_pos][1] = '\0';
								tokens[nonref_allele_column_pos][0] = 'R';
								tokens[nonref_allele_column_pos][1] = '\0';
							} else {
								log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
								if (!drop) write_columns();
								continue;
							}
						} else if (((auxiliary::strcmp_ignore_case(ref_allele, "D") != 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "R") != 0)) &&
									((auxiliary::strcmp_ignore_case(ref_allele, "R") != 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "D") != 0))) {
							log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
							if (!drop) write_columns();
							continue;
						}
					} else if (auxiliary::strcmp_ignore_case(buffer_tokens[2u], "I") == 0) {
						//CHR:POS:I
						query_position.type = 'I';
						if ((strspn(ref_allele, "ACGT") == ref_allele_length) && ((strspn(nonref_allele, "ACGT") == nonref_allele_length))) {
							if (ref_allele_length > nonref_allele_length) {
								log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "I", "R");
								tokens[ref_allele_column_pos][0] = 'I';
								tokens[ref_allele_column_pos][1] = '\0';
								tokens[nonref_allele_column_pos][0] = 'R';
								tokens[nonref_allele_column_pos][1] = '\0';
							} else if (ref_allele_length < nonref_allele_length) {
								log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "R", "I");
								tokens[ref_allele_column_pos][0] = 'R';
								tokens[ref_allele_column_pos][1] = '\0';
								tokens[nonref_allele_column_pos][0] = 'I';
								tokens[nonref_allele_column_pos][1] = '\0';
							} else {
								log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
								if (!drop) write_columns();
								continue;
							}
						} else if (((auxiliary::strcmp_ignore_case(ref_allele, "I") != 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "R") != 0)) &&
									((auxiliary::strcmp_ignore_case(ref_allele, "R") != 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "I") != 0))) {
							log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
							if (!drop) write_columns();
							continue;
						}
					}
				}

				if (query_position.type == '\0') {
					// CHR:POS:SOMETHING_ELSE or CHR:POS
					if (((auxiliary::strcmp_ignore_case(ref_allele, "A") == 0) || (auxiliary::strcmp_ignore_case(ref_allele, "C") == 0) ||
						(auxiliary::strcmp_ignore_case(ref_allele, "G") == 0) || (auxiliary::strcmp_ignore_case(ref_allele, "T") == 0)) &&
						((auxiliary::strcmp_ignore_case(nonref_allele, "A") == 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "C") == 0) ||
							(auxiliary::strcmp_ignore_case(nonref_allele, "G") == 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "T") == 0))) {
						query_position.type = 'S';
					} else if (((auxiliary::strcmp_ignore_case(ref_allele, "D") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "R") == 0)) ||
								((auxiliary::strcmp_ignore_case(ref_allele, "R") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "D") == 0))) {
						query_position.type = 'I';
					} else if (((auxiliary::strcmp_ignore_case(ref_allele, "I") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "R") == 0)) ||
								((auxiliary::strcmp_ignore_case(ref_allele, "R") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "I") == 0))) {
						query_position.type = 'I';
					} else {
						query_position.type = 'I';

						if (auxiliary::strcmp_ignore_case(ref_allele, nonref_allele) == 0) {
							log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
							if (!drop) write_columns();
							continue;
						}

						if ((ref_allele_length == nonref_allele_length) || (strspn(ref_allele, "ACGT") != ref_allele_length) || (strspn(nonref_allele, "ACGT") != nonref_allele_length)) {
							log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
							if (!drop) write_columns();
							continue;
						}

						recode_alleles = true;
					}
				}

				found_position_index_entry_pos = found_position - position_index;

				if (recode_alleles) {
					position = found_position_index_entry_pos;
					while ((position >= 0) && (qsort_position_index_entry_cmp(&(position_index[position]), &query_position) == 0)) {
						if (position_index[position].type == query_position.type) {
							consistent_type = true;

							if ((auxiliary::strcmp_ignore_case(position_index[position].ref_allele, ref_allele) == 0) &&
									(auxiliary::strcmp_ignore_case(position_index[position].nonref_allele, nonref_allele) == 0)) {
								consistent_alleles = true;

								if (ref_allele_length > nonref_allele_length) {
									log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "R", "D");
									tokens[ref_allele_column_pos][0] = 'R';
									tokens[ref_allele_column_pos][1] = '\0';
									tokens[nonref_allele_column_pos][0] = 'D';
									tokens[nonref_allele_column_pos][1] = '\0';
								} else if (ref_allele_length < nonref_allele_length) {
									log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "R", "I");
									tokens[ref_allele_column_pos][0] = 'R';
									tokens[ref_allele_column_pos][1] = '\0';
									tokens[nonref_allele_column_pos][0] = 'I';
									tokens[nonref_allele_column_pos][1] = '\0';
								}

								break;
							} else if ((auxiliary::strcmp_ignore_case(position_index[position].ref_allele, nonref_allele) == 0) &&
									(auxiliary::strcmp_ignore_case(position_index[position].nonref_allele, ref_allele) == 0)) {
								consistent_alleles = true;

								if (ref_allele_length > nonref_allele_length) {
									log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "I", "R");
									tokens[ref_allele_column_pos][0] = 'I';
									tokens[ref_allele_column_pos][1] = '\0';
									tokens[nonref_allele_column_pos][0] = 'R';
									tokens[nonref_allele_column_pos][1] = '\0';
								} else if (ref_allele_length < nonref_allele_length) {
									log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "D", "R");
									tokens[ref_allele_column_pos][0] = 'D';
									tokens[ref_allele_column_pos][1] = '\0';
									tokens[nonref_allele_column_pos][0] = 'R';
									tokens[nonref_allele_column_pos][1] = '\0';
								}

								break;
							}
						}
						--position;
					}

					if (!consistent_type) {
						position = found_position_index_entry_pos;
						++position;
						while ((position < index_size) && (qsort_position_index_entry_cmp(&(position_index[position]), &query_position) == 0)) {
							if (position_index[position].type == query_position.type) {
								consistent_type = true;

								if ((auxiliary::strcmp_ignore_case(position_index[position].ref_allele, ref_allele) == 0) &&
										(auxiliary::strcmp_ignore_case(position_index[position].nonref_allele, nonref_allele) == 0)) {
									consistent_alleles = true;

									if (ref_allele_length > nonref_allele_length) {
										log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "R", "D");
										tokens[ref_allele_column_pos][0] = 'R';
										tokens[ref_allele_column_pos][1] = '\0';
										tokens[nonref_allele_column_pos][0] = 'D';
										tokens[nonref_allele_column_pos][1] = '\0';
									} else if (ref_allele_length < nonref_allele_length) {
										log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "R", "I");
										tokens[ref_allele_column_pos][0] = 'R';
										tokens[ref_allele_column_pos][1] = '\0';
										tokens[nonref_allele_column_pos][0] = 'I';
										tokens[nonref_allele_column_pos][1] = '\0';
									}

									break;
								} else if ((auxiliary::strcmp_ignore_case(position_index[position].ref_allele, nonref_allele) == 0) &&
										(auxiliary::strcmp_ignore_case(position_index[position].nonref_allele, ref_allele) == 0)) {
									consistent_alleles = true;

									if (ref_allele_length > nonref_allele_length) {
										log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "I", "R");
										tokens[ref_allele_column_pos][0] = 'I';
										tokens[ref_allele_column_pos][1] = '\0';
										tokens[nonref_allele_column_pos][0] = 'R';
										tokens[nonref_allele_column_pos][1] = '\0';
									} else if (ref_allele_length < nonref_allele_length) {
										log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "D", "R");
										tokens[ref_allele_column_pos][0] = 'D';
										tokens[ref_allele_column_pos][1] = '\0';
										tokens[nonref_allele_column_pos][0] = 'R';
										tokens[nonref_allele_column_pos][1] = '\0';
									}

									break;
								}
							}
							++position;
						}
					}

					if (!consistent_type) {
						log_writer->write("Line %u: (WARNING) Type of %s doesn't match type in the map.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}

					if (!consistent_alleles) {
						log_writer->write("Line %u: (WARNING) Alleles of %s don't match alleles in the map.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}
				} else {
					position = found_position_index_entry_pos;
					while ((position >= 0) && (qsort_position_index_entry_cmp(&(position_index[position]), &query_position) == 0)) {
						if (position_index[position].type == query_position.type) {
							consistent_type = true;
							break;
						}
						--position;
					}

					if (!consistent_type) {
						position = found_position_index_entry_pos;
						++position;
						while ((position < index_size) && (qsort_position_index_entry_cmp(&(position_index[position]), &query_position) == 0)) {
							if (position_index[position].type == query_position.type) {
								consistent_type = true;
								break;
							}
							++position;
						}
					}

					if (!consistent_type) {
						log_writer->write("Line %u (WARNING): Type of %s doesn't match type in the map.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}
				}

				found_position = &(position_index[position]);
				if (found_position->type == 'S') {
					sprintf(buffer, "%s:%lu:%s", tokens[chr_column_pos], found_position->position, VCF_SNP_TYPE);
				} else {
					sprintf(buffer, "%s:%lu:%s", tokens[chr_column_pos], found_position->position, VCF_INDEL_TYPE_01);
				}

				if (strcmp(buffer,  tokens[id_column_pos]) != 0) {
					log_writer->write("Line %u: %s changed to %s.\n", line_number, tokens[id_column_pos], buffer);
					tokens[id_column_pos] = buffer;
				}
			}

			write_columns();
		}

		free(tokens);
		tokens = NULL;

		free(buffer);
		buffer = NULL;

		if (line_length == 0) {
			throw HarmonizerException("Harmonizer", "harmonize( bool )", __LINE__, 13, line_number + 1u, input_file);
		}
	} catch (ReaderException &e) {
		HarmonizerException new_e(e);
		new_e.add_message("Harmonizer", "harmonize( bool )", __LINE__, 5, (input_file != NULL) ? input_file : "NULL");
		throw new_e;
	} catch (WriterException &e) {
		HarmonizerException new_e(e);
		new_e.add_message("Harmonizer", "harmonize( bool )", __LINE__, 5, (input_file != NULL) ? input_file : "NULL");
		throw new_e;
	} catch (HarmonizerException &e) {
		e.add_message("Harmonizer", "harmonize( bool )", __LINE__, 5, (input_file != NULL) ? input_file : "NULL");
		throw;
	}
}

void Harmonizer::harmonize_no_vcf_allele_check(bool drop) throw (HarmonizerException) {
	char* line = NULL;
	int line_length = 0;
	unsigned int line_number = 1u;

	char* token = NULL;

	int column_number = 0;

	chr_index* index = NULL;
	unsigned int index_size = 0u;

	position_index_entry* position_index = NULL;
	position_index_entry query_position;
	position_index_entry* found_position = NULL;

	id_index_entry* id_index = NULL;
	id_index_entry query_id;
	id_index_entry* found_id = NULL;

	unsigned int found_id_index_entry_pos = 0u;
	unsigned int found_position_index_entry_pos = 0u;
	long int position = 0;

	char* ref_allele = NULL;
	char* nonref_allele = NULL;
	unsigned int ref_allele_length = 0u;
	unsigned int nonref_allele_length = 0u;

	char* buffer = NULL;
	char* buffer_tmp = NULL;
	unsigned int buffer_tokens_number = 0u;
	char* buffer_token = NULL;
	char* buffer_tokens[3u];
	char* end_ptr = NULL;

	bool empty_ref_allele = false;
	bool empty_nonref_allele = false;
	bool consistent_type = false;

	if ((input_file == NULL) || (reader == NULL) ||
			(output_file == NULL) || (writer == NULL) ||
			(log_file == NULL) || (log_writer == NULL) ||
			(chr_column == NULL) || (chr_column_pos < 0) ||
			(id_column == NULL) || (id_column_pos < 0) ||
			(ref_allele_column == NULL) || (ref_allele_column_pos < 0) ||
			(nonref_allele_column == NULL) || (nonref_allele_column_pos < 0)) {
		return;
	}

	try {
		tokens = (char**)malloc(file_column_number * sizeof(char*));
		if (tokens == NULL) {
			throw HarmonizerException("Harmonizer", "harmonize( bool )", __LINE__, 2, file_column_number * sizeof(char*));
		}

		buffer = (char*)malloc(16777216 * sizeof(char));
		if (buffer == NULL) {
			throw HarmonizerException("Harmonizer", "harmonize( bool )", __LINE__, 2, 16777216 * sizeof(char));
		}

		writer->write("%s\n", header_backup);

		while ((line_length = reader->read_line()) > 0) {
			line = *(reader->line);
			++line_number;

			consistent_type = false;

			column_number = 0;
			while ((token = auxiliary::strtok(&line, separator)) != NULL) {
				if (column_number < file_column_number) {
					tokens[column_number] = token;
				}
				++column_number;
			}

			if (column_number < file_column_number) {
				throw HarmonizerException("Harmonizer", "harmonize( bool )", __LINE__, 14, line_number, input_file, column_number, file_column_number);
			} else if (column_number > file_column_number) {
				throw HarmonizerException("Harmonizer", "harmonize( bool )", __LINE__, 15, line_number, input_file, column_number, file_column_number);
			}

			auxiliary::trim(&(tokens[chr_column_pos]));
			map_index_by_chr_it = map_index_by_chr.find(tokens[chr_column_pos]);
			if (map_index_by_chr_it == map_index_by_chr.end()) {
				log_writer->write("Line %u: (WARNING) Chromosome %s of %s is not in the map.\n", line_number, tokens[chr_column_pos], tokens[id_column_pos]);
				if (!drop) write_columns();
				continue;
			}

			index = map_index_by_chr_it->second;
			index_size = index->n;

			position_index = index->positions;
			id_index = index->ids;

			found_id = NULL;
			found_position = NULL;

			auxiliary::trim(&(tokens[ref_allele_column_pos]));
			ref_allele = tokens[ref_allele_column_pos];
			ref_allele_length = strlen(ref_allele);

			auxiliary::trim(&(tokens[nonref_allele_column_pos]));
			nonref_allele = tokens[nonref_allele_column_pos];
			nonref_allele_length = strlen(nonref_allele);

			auxiliary::trim(&(tokens[id_column_pos]));
			if (auxiliary::strcmp_ignore_case(tokens[id_column_pos], "rs", 2) == 0) { // Look up by rsID
				query_id.id = tokens[id_column_pos];
				found_id = (id_index_entry*)bsearch(&query_id, id_index, index_size, sizeof(id_index_entry), qsort_id_index_entry_cmp);
				if (found_id == NULL) {
					log_writer->write("Line %u: (WARNING) %s is not in the map.\n", line_number, query_id.id);
					if (!drop) write_columns();
					continue;
				}

				if (((auxiliary::strcmp_ignore_case(ref_allele, "A") == 0) || (auxiliary::strcmp_ignore_case(ref_allele, "C") == 0) ||
					(auxiliary::strcmp_ignore_case(ref_allele, "G") == 0) || (auxiliary::strcmp_ignore_case(ref_allele, "T") == 0)) &&
					((auxiliary::strcmp_ignore_case(nonref_allele, "A") == 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "C") == 0) ||
						(auxiliary::strcmp_ignore_case(nonref_allele, "G") == 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "T") == 0))) {
					query_id.type = 'S';
				} else if (((auxiliary::strcmp_ignore_case(ref_allele, "D") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "R") == 0)) ||
							((auxiliary::strcmp_ignore_case(ref_allele, "R") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "D") == 0))) {
					query_id.type = 'I';
				} else if (((auxiliary::strcmp_ignore_case(ref_allele, "I") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "R") == 0)) ||
							((auxiliary::strcmp_ignore_case(ref_allele, "R") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "I") == 0))) {
					query_id.type = 'I';
				} else {
					query_id.type = 'I';

					if (auxiliary::strcmp_ignore_case(ref_allele, nonref_allele) == 0) {
						log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}

					if ((ref_allele_length == nonref_allele_length) || (strspn(ref_allele, "ACGT") != ref_allele_length) || (strspn(nonref_allele, "ACGT") != nonref_allele_length)) {
						log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}
				}

				found_id_index_entry_pos = found_id - id_index;

				position = found_id_index_entry_pos;
				while ((position >= 0) && (qsort_id_index_entry_cmp(&(id_index[position]), &query_id) == 0)) {
					if (position_index[id_index[position].location].type == query_id.type) {
						consistent_type = true;
						break;
					}
					--position;
				}

				if (!consistent_type) {
					position = found_id_index_entry_pos;
					++position;
					while ((position < index_size) && (qsort_id_index_entry_cmp(&(id_index[position]), &query_id) == 0)) {
						if (position_index[id_index[position].location].type == query_id.type) {
							consistent_type = true;
							break;
						}
						++position;
					}
				}

				if (!consistent_type) {
					log_writer->write("Line %u: (WARNING) Type of %s doesn't match type in the map.\n", line_number, query_id.id);
					if (!drop) write_columns();
					continue;
				}

				found_position = &(position_index[id_index[position].location]);
				if (found_position->type == 'S') {
					sprintf(buffer, "%s:%lu:%s", tokens[chr_column_pos], found_position->position, VCF_SNP_TYPE);
				} else {
					sprintf(buffer, "%s:%lu:%s", tokens[chr_column_pos], found_position->position, VCF_INDEL_TYPE_01);
				}
				log_writer->write("Line %u: %s changed to %s.\n", line_number, tokens[id_column_pos], buffer);
				tokens[id_column_pos] = buffer;
			} else if (auxiliary::strcmp_ignore_case(tokens[id_column_pos], "merged_del", 10) == 0) { // Look up by id
				query_id.id = tokens[id_column_pos];
				found_id = (id_index_entry*)bsearch(&query_id, id_index, index_size, sizeof(id_index_entry), qsort_id_index_entry_cmp);
				if (found_id == NULL) {
					log_writer->write("Line %u: (WARNING) %s is not in the map.\n", line_number, query_id.id);
					if (!drop) write_columns();
					continue;
				}

				if (((auxiliary::strcmp_ignore_case(ref_allele, "D") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "R") == 0)) ||
							((auxiliary::strcmp_ignore_case(ref_allele, "R") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "D") == 0))) {
					query_id.type = 'I';
				} else {
					query_id.type = 'I';

					empty_ref_allele = false;
					empty_nonref_allele = false;

					if ((strcmp(ref_allele, "-") == 0) || (strcmp(ref_allele, ".") == 0) || (auxiliary::strcmp_ignore_case(ref_allele, "NA") == 0)) {
						empty_ref_allele = true;
					} else if (strspn(ref_allele, "ACGT") != ref_allele_length) {
						log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}

					if ((strcmp(nonref_allele, "-") == 0) || (strcmp(nonref_allele, ".") == 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "NA") == 0)) {
						empty_nonref_allele = true;
					} else if (strspn(nonref_allele, "ACGT") != nonref_allele_length) {
						log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}

					if (empty_ref_allele && empty_nonref_allele) {
						log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}

					if (!empty_ref_allele && !empty_nonref_allele && (ref_allele_length == nonref_allele_length)) {
						log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}

					if (empty_ref_allele) {
						log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "D", "R");
						tokens[ref_allele_column_pos][0] = 'D';
						tokens[ref_allele_column_pos][1] = '\0';
						tokens[nonref_allele_column_pos][0] = 'R';
						tokens[nonref_allele_column_pos][1] = '\0';
					} else if (empty_nonref_allele) {
						log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "R", "D");
						tokens[ref_allele_column_pos][0] = 'R';
						tokens[ref_allele_column_pos][1] = '\0';
						tokens[nonref_allele_column_pos][0] = 'D';
						tokens[nonref_allele_column_pos][1] = '\0';
					} else if (ref_allele_length > nonref_allele_length) {
						log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "R", "D");
						tokens[ref_allele_column_pos][0] = 'R';
						tokens[ref_allele_column_pos][1] = '\0';
						tokens[nonref_allele_column_pos][0] = 'D';
						tokens[nonref_allele_column_pos][1] = '\0';
					} else if (ref_allele_length < nonref_allele_length) {
						log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "D", "R");
						tokens[ref_allele_column_pos][0] = 'D';
						tokens[ref_allele_column_pos][1] = '\0';
						tokens[nonref_allele_column_pos][0] = 'R';
						tokens[nonref_allele_column_pos][1] = '\0';
					} else {
						log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
						if (!drop) write_columns();
						continue;
					}
				}

				found_id_index_entry_pos = found_id - id_index;

				position = found_id_index_entry_pos;
				while ((position >= 0) && (qsort_id_index_entry_cmp(&(id_index[position]), &query_id) == 0)) {
					if (position_index[id_index[position].location].type == query_id.type) {
						consistent_type = true;
						break;
					}
					--position;
				}

				if (!consistent_type) {
					position = found_id_index_entry_pos;
					++position;
					while ((position < index_size) && (qsort_id_index_entry_cmp(&(id_index[position]), &query_id) == 0)) {
						if (position_index[id_index[position].location].type == query_id.type) {
							consistent_type = true;
							break;
						}
						++position;
					}
				}

				if (!consistent_type) {
					log_writer->write("Line %u: (WARNING) Type of %s doesn't match type in the map.\n", line_number, query_id.id);
					if (!drop) write_columns();
					continue;
				}

				found_position = &(position_index[id_index[position].location]);
				if (found_position->type == 'S') {
					sprintf(buffer, "%s:%lu:%s", tokens[chr_column_pos], found_position->position, VCF_SNP_TYPE);
				} else {
					sprintf(buffer, "%s:%lu:%s", tokens[chr_column_pos], found_position->position, VCF_INDEL_TYPE_01);
				}
				log_writer->write("Line %u: %s changed to %s.\n", line_number, tokens[id_column_pos], buffer);
				tokens[id_column_pos] = buffer;
			} else { // Look up by position
				strcpy(buffer, tokens[id_column_pos]);
				buffer_tmp = buffer;

				buffer_tokens_number = 0u;
				buffer_tokens[0u] = buffer_tokens[1u] = buffer_tokens[2u] = NULL;
				while ((buffer_tokens_number < 3u) && ((buffer_token = auxiliary::strtok(&buffer_tmp, ':')) != NULL)) {
					buffer_tokens[buffer_tokens_number++] = buffer_token;
				}

				if ((buffer_tokens_number < 2u) || (buffer_tokens_number > 3u)) {
					log_writer->write("Line %u: (WARNING) %s identifier doesn't follow CHR:POS, CHR:POS:TYPE, or rsID formats.\n", line_number, tokens[id_column_pos]);
					if (!drop) write_columns();
					continue;
				}

				query_position.position = strtoul(buffer_tokens[1u], &end_ptr, 10);
				if (end_ptr == buffer_tokens[1u]) {
					log_writer->write("Line %u: (WARNING) Position of %s can't be parsed to integer.\n", line_number, tokens[id_column_pos]);
					if (!drop) write_columns();
					continue;
				}

				found_position = (position_index_entry*)bsearch(&query_position, position_index, index_size, sizeof(position_index_entry), qsort_position_index_entry_cmp);
				if (found_position == NULL) {
					log_writer->write("Line %u: (WARNING) %s is not in the map.\n", line_number, tokens[id_column_pos]);
					if (!drop) write_columns();
					continue;
				}

				query_position.type = '\0';
				if (buffer_tokens[2u] != NULL) {
					if (auxiliary::strcmp_ignore_case(buffer_tokens[2u], VCF_SNP_TYPE) == 0) {
						//CHR:POS:SNP
						query_position.type = 'S';
						if (((auxiliary::strcmp_ignore_case(ref_allele, "A") != 0) && (auxiliary::strcmp_ignore_case(ref_allele, "C") != 0) &&
								(auxiliary::strcmp_ignore_case(ref_allele, "G") != 0) && (auxiliary::strcmp_ignore_case(ref_allele, "T") != 0)) ||
								((auxiliary::strcmp_ignore_case(nonref_allele, "A") != 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "C") != 0) &&
										(auxiliary::strcmp_ignore_case(nonref_allele, "G") != 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "T") != 0))) {
							log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
							if (!drop) write_columns();
							continue;
						}
					} else if (auxiliary::strcmp_ignore_case(buffer_tokens[2u], "D") == 0) {
						//CHR:POS:D
						query_position.type = 'I';
						if ((strspn(ref_allele, "ACGT") == ref_allele_length) && ((strspn(nonref_allele, "ACGT") == nonref_allele_length))) {
							if (ref_allele_length > nonref_allele_length) {
								log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "R", "D");
								tokens[ref_allele_column_pos][0] = 'R';
								tokens[ref_allele_column_pos][1] = '\0';
								tokens[nonref_allele_column_pos][0] = 'D';
								tokens[nonref_allele_column_pos][1] = '\0';
							} else if (ref_allele_length < nonref_allele_length) {
								log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "D", "R");
								tokens[ref_allele_column_pos][0] = 'D';
								tokens[ref_allele_column_pos][1] = '\0';
								tokens[nonref_allele_column_pos][0] = 'R';
								tokens[nonref_allele_column_pos][1] = '\0';
							} else {
								log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
								if (!drop) write_columns();
								continue;
							}
						} else if (((auxiliary::strcmp_ignore_case(ref_allele, "D") != 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "R") != 0)) &&
									((auxiliary::strcmp_ignore_case(ref_allele, "R") != 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "D") != 0))) {
							log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
							if (!drop) write_columns();
							continue;
						}
					} else if (auxiliary::strcmp_ignore_case(buffer_tokens[2u], "I") == 0) {
						//CHR:POS:I
						query_position.type = 'I';
						if ((strspn(ref_allele, "ACGT") == ref_allele_length) && ((strspn(nonref_allele, "ACGT") == nonref_allele_length))) {
							if (ref_allele_length > nonref_allele_length) {
								log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "I", "R");
								tokens[ref_allele_column_pos][0] = 'I';
								tokens[ref_allele_column_pos][1] = '\0';
								tokens[nonref_allele_column_pos][0] = 'R';
								tokens[nonref_allele_column_pos][1] = '\0';
							} else if (ref_allele_length < nonref_allele_length) {
								log_writer->write("Line %u: Alleles %s/%s of %s changed to %s/%s.\n", line_number, ref_allele, nonref_allele, tokens[id_column_pos], "R", "I");
								tokens[ref_allele_column_pos][0] = 'R';
								tokens[ref_allele_column_pos][1] = '\0';
								tokens[nonref_allele_column_pos][0] = 'I';
								tokens[nonref_allele_column_pos][1] = '\0';
							} else {
								log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
								if (!drop) write_columns();
								continue;
							}
						} else if (((auxiliary::strcmp_ignore_case(ref_allele, "I") != 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "R") != 0)) &&
									((auxiliary::strcmp_ignore_case(ref_allele, "R") != 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "I") != 0))) {
							log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
							if (!drop) write_columns();
							continue;
						}
					}
				}

				if (query_position.type == '\0') {
					// CHR:POS:SOMETHING_ELSE or CHR:POS
					if (((auxiliary::strcmp_ignore_case(ref_allele, "A") == 0) || (auxiliary::strcmp_ignore_case(ref_allele, "C") == 0) ||
						(auxiliary::strcmp_ignore_case(ref_allele, "G") == 0) || (auxiliary::strcmp_ignore_case(ref_allele, "T") == 0)) &&
						((auxiliary::strcmp_ignore_case(nonref_allele, "A") == 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "C") == 0) ||
							(auxiliary::strcmp_ignore_case(nonref_allele, "G") == 0) || (auxiliary::strcmp_ignore_case(nonref_allele, "T") == 0))) {
						query_position.type = 'S';
					} else if (((auxiliary::strcmp_ignore_case(ref_allele, "D") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "R") == 0)) ||
								((auxiliary::strcmp_ignore_case(ref_allele, "R") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "D") == 0))) {
						query_position.type = 'I';
					} else if (((auxiliary::strcmp_ignore_case(ref_allele, "I") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "R") == 0)) ||
								((auxiliary::strcmp_ignore_case(ref_allele, "R") == 0) && (auxiliary::strcmp_ignore_case(nonref_allele, "I") == 0))) {
						query_position.type = 'I';
					} else {
						query_position.type = 'I';

						if (auxiliary::strcmp_ignore_case(ref_allele, nonref_allele) == 0) {
							log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
							if (!drop) write_columns();
							continue;
						}

						if ((ref_allele_length == nonref_allele_length) || (strspn(ref_allele, "ACGT") != ref_allele_length) || (strspn(nonref_allele, "ACGT") != nonref_allele_length)) {
							log_writer->write("Line %u: (WARNING) Alleles of %s don't match the specified variation type.\n", line_number, tokens[id_column_pos]);
							if (!drop) write_columns();
							continue;
						}
					}
				}

				found_position_index_entry_pos = found_position - position_index;

				position = found_position_index_entry_pos;
				while ((position >= 0) && (qsort_position_index_entry_cmp(&(position_index[position]), &query_position) == 0)) {
					if (position_index[position].type == query_position.type) {
						consistent_type = true;
						break;
					}
					--position;
				}

				if (!consistent_type) {
					position = found_position_index_entry_pos;
					++position;
					while ((position < index_size) && (qsort_position_index_entry_cmp(&(position_index[position]), &query_position) == 0)) {
						if (position_index[position].type == query_position.type) {
							consistent_type = true;
							break;
						}
						++position;
					}
				}

				if (!consistent_type) {
					log_writer->write("Line %u (WARNING): Type of %s doesn't match type in the map.\n", line_number, tokens[id_column_pos]);
					if (!drop) write_columns();
					continue;
				}

				found_position = &(position_index[position]);
				if (found_position->type == 'S') {
					sprintf(buffer, "%s:%lu:%s", tokens[chr_column_pos], found_position->position, VCF_SNP_TYPE);
				} else {
					sprintf(buffer, "%s:%lu:%s", tokens[chr_column_pos], found_position->position, VCF_INDEL_TYPE_01);
				}

				if (strcmp(buffer,  tokens[id_column_pos]) != 0) {
					log_writer->write("Line %u: %s changed to %s.\n", line_number, tokens[id_column_pos], buffer);
					tokens[id_column_pos] = buffer;
				}
			}

			write_columns();
		}

		free(tokens);
		tokens = NULL;

		free(buffer);
		buffer = NULL;

		if (line_length == 0) {
			throw HarmonizerException("Harmonizer", "harmonize( bool )", __LINE__, 13, line_number + 1u, input_file);
		}
	} catch (ReaderException &e) {
		HarmonizerException new_e(e);
		new_e.add_message("Harmonizer", "harmonize( bool )", __LINE__, 5, (input_file != NULL) ? input_file : "NULL");
		throw new_e;
	} catch (WriterException &e) {
		HarmonizerException new_e(e);
		new_e.add_message("Harmonizer", "harmonize( bool )", __LINE__, 5, (input_file != NULL) ? input_file : "NULL");
		throw new_e;
	} catch (HarmonizerException &e) {
		e.add_message("Harmonizer", "harmonize( bool )", __LINE__, 5, (input_file != NULL) ? input_file : "NULL");
		throw;
	}
}

inline int Harmonizer::qsort_position_index_entry_cmp(const void* first, const void* second) {
	position_index_entry* first_entry = (position_index_entry*)first;
	position_index_entry* second_entry = (position_index_entry*)second;

	if (first_entry->position > second_entry->position) {
		return 1;
	} else if (first_entry->position < second_entry->position) {
		return -1;
	}

	return 0;
}

inline int Harmonizer::qsort_id_index_entry_cmp(const void* first, const void* second) {
	return auxiliary::strcmp_ignore_case(((id_index_entry*)first)->id, ((id_index_entry*)second)->id);
}
