/**************************************/
#pragma once
/**************************************/
#include <stdio.h>
#include <stdint.h>
/**************************************/

//! FourCC from string
#define RIFF_FOURCC(x)	( \
			((uint32_t)(x)[0]&0xFF)       | \
			((uint32_t)(x)[1]&0xFF) <<  8 | \
			((uint32_t)(x)[2]&0xFF) << 16 | \
			((uint32_t)(x)[3]&0xFF) << 24   \
			)

//! RIFF chunk header structure
struct RIFF_CkHeader_t {
	uint32_t Type;
	uint32_t Size;
};

/**************************************/

//! RIFF chunk handler structure
//! Finish a list of these by appending (struct RIFF_CkHdl_t){.Type=0}
struct RIFF_CkHdl_t {
	uint32_t Type;
	int (*Func)(FILE *File, void *User, const struct RIFF_CkHeader_t *Ck);
};

//! RIFF 'LIST' chunk handler structure
//! Finish a list of these by appending (struct RIFF_CkListHdl_t){.Type=0}
struct RIFF_CkListHdl_t {
	uint32_t Type;
	const struct RIFF_CkHdl_t     *CkHdl;
	const struct RIFF_CkListHdl_t *ListHdl;
	int (*ListCbBeg)(FILE *File, void *User);
	int (*ListCbEnd)(FILE *File, void *User);
};

/**************************************/

//! RIFF_CkRead(File, User, CkHdl, ListHdl)
//! Description: Read RIFF chunk.
//! Arguments:
//!   File:    Source file.
//!   User:    Userdata to pass to handlers.
//!   CkHdl:   List of chunk handlers.
//!   ListHdl: List of RIFF/LIST handlers.
//! Returns:
//!   On no handler found, returns 0.
//!   Otherwise, returns the value of the last-handled chunk.
//!   When encounternig a RIFF/LIST chunk, execution stops when
//!   an associated handler returns a value less than 0.
int RIFF_CkRead(FILE *File, void *User, const struct RIFF_CkHdl_t *CkHdl, const struct RIFF_CkListHdl_t *ListHdl);

/**************************************/
//! EOF
/**************************************/
