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

#include "include/Harmonizer2Exception.h"

const int Harmonizer2Exception::MESSAGE_TEMPLATES_NUMBER = 22;
const char* Harmonizer2Exception::MESSAGE_TEMPLATES[] = {
/*00*/	"The '%s' argument has NULL value.",
/*01*/	"The '%s' argument has an invalid value.",
/*02*/	"Memory allocation error (%d bytes).",
/*03*/	"Error while opening '%s' file.",
/*04*/	"Error while closing '%s' file.",
/*05*/	"Error while reading '%s' file.",
/*06*/	"Error while writing '%s' file.",
/*07*/	"Error while indexing VCF file.",
/*08*/	"The mandatory VCF file format information line is incorrect.",
/*09*/	"The mandatory VCF file format information line is empty.",
/*10*/	"The mandatory VCF header line was not found.",
/*11*/	"The VCF file column '%s' is missing on position %d.",
/*12*/	"The VCF header/meta-information on line %d is incorrect.",
/*13*/	"VCF file line %d is empty.",
/*14*/	"Line %d in '%s' file contains less columns (%d) than expected (%d).",
/*15*/	"Line %d in '%s' file contains more columns (%d) than expected (%d).",
/*16*/	"The header is missing on line %d in '%s' file.",
/*17*/	"Column '%s' is missing in the header of '%s' file.",
/*18*/	"Error while initializing a log file name.",
/*19*/	"Memory re-allocation error (%d bytes).",
/*20*/	"The chromosomal position '%s' on line %d in file '%s' could not be parsed to unsigned integer."
/*21*/	"Line %d in '%s' file is empty."
};

Harmonizer2Exception::Harmonizer2Exception() : Exception() {

}

Harmonizer2Exception::Harmonizer2Exception(int message_template_index, ... ) : Exception()  {
	va_list arguments;

	va_start(arguments, message_template_index);
	add_message(message_template_index, arguments);
	va_end(arguments);
}

Harmonizer2Exception::Harmonizer2Exception(const char* class_name, const char* method_name, int source_line) : Exception(class_name, method_name, source_line)  {

}

Harmonizer2Exception::Harmonizer2Exception(const char* class_name, const char* method_name, int source_line, int message_template_index, ... ) : Exception() {
	va_list arguments;

	va_start(arguments, message_template_index);
	add_message(class_name, method_name, source_line, message_template_index, arguments);
	va_end(arguments);
}

Harmonizer2Exception::Harmonizer2Exception(const Exception& exception) : Exception(exception) {

}

Harmonizer2Exception::~Harmonizer2Exception() throw() {

}

const char* Harmonizer2Exception::get_message_template(int message_template_index) {
	if ((message_template_index >= 0) && (message_template_index < MESSAGE_TEMPLATES_NUMBER)) {
		return  MESSAGE_TEMPLATES[message_template_index];
	}

	return NULL;
}
