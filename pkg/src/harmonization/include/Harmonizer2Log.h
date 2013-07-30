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

#ifndef HARMONIZER2LOG_H_
#define HARMONIZER2LOG_H_

#include "Harmonizer2Exception.h"

class Harmonizer2Log {
public:
	enum warning_type {
		CHROMOSOME_MISSING,
		POSITION_MISSING,
		TYPE_MISMATCH,
		ALLELE_MISMATCH,
		STRAND_MISMATCH,
		BAD_POSITION,
		BAD_ALLELES
	};

	enum message_type {
		STRAND_FLIPPED,
		ALLELES_CHANGED,
		ID_CHANGED
	};

	static const int N_WARNING_TYPES;
	static const char* warnings[];
	static const int N_MESSAGE_TYPES;
	static const char* messages[];

	unsigned int* warning_counts;
	unsigned int* message_counts;

	Harmonizer2Log() throw (Harmonizer2Exception);
	virtual ~Harmonizer2Log();

	void add_warning(warning_type type);
	void add_message(message_type type);
};

#endif
