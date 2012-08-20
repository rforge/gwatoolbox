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

#include "include/GzipReader.h"

const unsigned int GzipReader::DEFAULT_BUFFER_SIZE = 16777216;

GzipReader::GzipReader(unsigned int buffer_size) throw (ReaderException) : Reader(&buffer),
	opened(false), buffer_size(buffer_size), buffer(NULL) {

	buffer = (char*)malloc((buffer_size + 1) * sizeof(char));
	if (buffer == NULL) {
		throw ReaderException("GzipReader", "GzipReader( int )", __LINE__, 2, (buffer_size + 1) * sizeof(char));
	}

	buffer[0] = '\0';
}

GzipReader::~GzipReader() {
	buffer_size = 0;

	free(buffer);
	buffer = NULL;
}

void GzipReader::open() throw (ReaderException) {
	close();

	infile = gzopen(file_name, "rb");
	if (infile == NULL) {
		throw ReaderException("GzipReader", "open()", __LINE__, 3, file_name);
	}

	opened = true;
}

void GzipReader::close() throw (ReaderException) {
	if (opened) {
		int gzerrno = 0;

		gzerrno = gzclose(infile);
		if (gzerrno != Z_OK) {
			throw ReaderException("GzipReader", "open()", __LINE__, 5, file_name);
		}

		opened = false;
	}
}

int GzipReader::read_line() throw (ReaderException) {
	int i = 0;
	int c = 0;

	while ((i < buffer_size) && ((c = gzgetc(infile)) >= 0)) {
		buffer[i] = (char)c;

		if (buffer[i] == '\n') {
			buffer[i] = '\0';
			return i;
		} else if (buffer[i] == '\r') {
			buffer[i] = '\0';
			if ((c = gzgetc(infile)) >= 0) {
				if ((char)c != '\n') {
					c = gzungetc(c, infile);
				}
			}
			return i;
		}

		i += 1;
	}

	buffer[i] = '\0';

	if ((c < 0) && (gzeof(infile) < 1)) {
		throw ReaderException("GzipReader", "int read_line()", __LINE__, 4, file_name);
	}

	return (i == 0 ? -1 : i);
}

void GzipReader::reset() throw (ReaderException) {
	if (gzseek(infile, 0L, SEEK_SET) < 0) {
		throw ReaderException("GzipReader", "reset()", __LINE__, 6, file_name);
	}
}

bool GzipReader::eof() {
	return gzeof(infile) > 0;
}

bool GzipReader::is_open() {
	return opened;
}

bool GzipReader::is_compressed() {
	return true;
}
