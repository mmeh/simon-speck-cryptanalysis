// ===========================================================================
// SIMON implementation and cryptanalytic methods
// =========================================================================
// Copyright (c) 2013 Martin M. Lauridsen and Hoda A. Alkhzaimi.

// Permission to use, copy, modify, and/or distribute this software for any
// purpose with or without fee is hereby granted, provided that the above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
// WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
// ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
// ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
// OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

#ifndef _SIMON_H_
#define _SIMON_H_

#include <stdint.h>
#include <string>

// Cipher setup parameters
#define BLOCK_SIZE (16)
#define KEY_WORDS (4)
#define PRINT_ROUND_KEYS (0)

// Rest of parameters are automated
#if (BLOCK_SIZE == 64)
	#define WORD_MASK (0xffffffffffffffffull)
#else
	#define WORD_MASK ((0x1ull << (BLOCK_SIZE&63)) - 1)
#endif
#define CONST_C ((0xffffffffffffffffull ^ 0x3ull) & WORD_MASK)

#if (BLOCK_SIZE == 4)
	#define ROUNDS (32)
	#define CONST_J (0)
#elif (BLOCK_SIZE == 16)
	#define ROUNDS (32)
	#define CONST_J (0)
#elif (BLOCK_SIZE == 24)
	#if (KEY_WORDS == 3)
		#define ROUNDS (36)
		#define CONST_J (0)
	#elif (KEY_WORDS == 4)
		#define ROUNDS (36)
		#define CONST_J (1)
	#endif
#elif (BLOCK_SIZE == 32)
	#if (KEY_WORDS == 3)
		#define ROUNDS (42)
		#define CONST_J (2)
	#elif (KEY_WORDS == 4)
		#define ROUNDS (44)
		#define CONST_J (3)
	#endif
#elif (BLOCK_SIZE == 48)
	#if (KEY_WORDS == 2)
		#define ROUNDS (52)
		#define CONST_J (2)
	#elif (KEY_WORDS == 3)
		#define ROUNDS (54)
		#define CONST_J (3)
	#endif
#elif (BLOCK_SIZE == 64)
	#if (KEY_WORDS == 2)
		#define ROUNDS (68)
		#define CONST_J (2)
	#elif (KEY_WORDS == 3)
		#define ROUNDS (69)
		#define CONST_J (3)
	#elif (KEY_WORDS == 4)
		#define ROUNDS (72)
		#define CONST_J (4)
	#endif
#endif

#define STATE_MASK ((WORD_MASK << BLOCK_SIZE) | WORD_MASK)

typedef uint64_t u64;

extern u64 z[5][62];
extern u64 k[ROUNDS];

// HELPERS
char* binary(u64);
int weight(u64);
unsigned char lcs(char *, char *);

// CIPHER STUFF
u64 rotate(u64, int);
u64 F(u64);
u64 F_toy(u64);
void key_schedule();
void encrypt(u64 &, u64 &, int);
void encrypt(u64 &, u64 &);
void decrypt(u64 &, u64 &, int);
void decrypt(u64 &, u64 &);
void setup_random_key();

// TESTING
void run_test_vectors();
void test_enc();

// DIFFERENTIAL STUFF
void distribution_test();
void diff_dist_table();
void test_differential(u64, u64, u64, u64, int);
void diff_attack(u64, u64, u64, u64, int);
void test_2r_difference(u64);
void ddt_diagonal();
void differences_to_zero_diff();

// ROTATIONAL STUFF
void rotational_approx();

// WEAK KEYS STUFF
void weak_keys();
void test_rotational();
void test_key_rotation();
void is_F_balanced();
void generate_key_relations(unsigned int);
void key_difference();

// IMPOSSIBLE DIFFERENTIAL STUFF
std::string rotchar(std::string, unsigned int);
unsigned short cweight(std::string);
std::string cxor(std::string, std::string);
std::string rot_func(std::string);
void impossible_diff_attack();
int impossible_diff(std::string, std::string);
void imp_diff_attack2();

#endif