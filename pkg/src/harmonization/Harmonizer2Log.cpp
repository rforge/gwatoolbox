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

#include "include/Harmonizer2Log.h"

const int Harmonizer2Log::N_WARNING_TYPES = 7;
const char* Harmonizer2Log::warnings[N_WARNING_TYPES] = {
		"Chromosome name is not in VCF",
		"Position is not in VCF",
		"Variation type mismatch in input file and VCF",
		"Alleles mismatch in input file and VCF",
		"Strand mismatch in input file and VCF",
		"Position is not numeric in input file",
		"Unrecognized alleles in input file"
};
const int Harmonizer2Log::N_MESSAGE_TYPES = 3;
const char* Harmonizer2Log::messages[N_MESSAGE_TYPES] = {
		"Strand flipped",
		"Alleles changed",
		"ID changed"
};

Harmonizer2Log::Harmonizer2Log() throw (Harmonizer2Exception) : warning_counts(NULL), message_counts(NULL) {
	warning_counts = (unsigned int*)malloc(N_WARNING_TYPES * sizeof(unsigned int));
	if (warning_counts == NULL) {
		throw Harmonizer2Exception("Harmonizer2Log", "Harmonizer2Log()", __LINE__, 2, (N_WARNING_TYPES * sizeof(unsigned int)));
	}

	message_counts = (unsigned int*)malloc(N_MESSAGE_TYPES * sizeof(unsigned int));
	if (message_counts == NULL) {
		throw Harmonizer2Exception("Harmonizer2Log", "Harmonizer2Log()", __LINE__, 2, (N_MESSAGE_TYPES * sizeof(unsigned int)));
	}

	for (int i = 0; i < N_WARNING_TYPES; ++i) {
		warning_counts[i] = 0;
	}

	for (int i = 0; i < N_MESSAGE_TYPES; ++i) {
		message_counts[i] = 0;
	}
}

Harmonizer2Log::~Harmonizer2Log() {
	if (warning_counts != NULL) {
		free(warning_counts);
		warning_counts = NULL;
	}

	if (message_counts != NULL) {
		free(message_counts);
		message_counts = NULL;
	}
}

void Harmonizer2Log::add_warning(warning_type type) {
	warning_counts[type] += 1u;
}

void Harmonizer2Log::add_message(message_type type) {
	message_counts[type] += 1u;
}
