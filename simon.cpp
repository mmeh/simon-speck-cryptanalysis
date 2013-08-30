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

#include "simon.h"
#include <stdio.h>
#include <string.h>
#include <random>
#include <vector>
#include <map>
#include <set>
#include "omp.h"

u64 z[5][62] = {
	{1,1,1,1,1,0,1,0,0,0,1,0,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,0,1,1,0,1,1,1,1,1,0,1,0,0,0,1,0,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,0,1,1,0},
	{1,0,0,0,1,1,1,0,1,1,1,1,1,0,0,1,0,0,1,1,0,0,0,0,1,0,1,1,0,1,0,1,0,0,0,1,1,1,0,1,1,1,1,1,0,0,1,0,0,1,1,0,0,0,0,1,0,1,1,0,1,0},
	{1,0,1,0,1,1,1,1,0,1,1,1,0,0,0,0,0,0,1,1,0,1,0,0,1,0,0,1,1,0,0,0,1,0,1,0,0,0,0,1,0,0,0,1,1,1,1,1,1,0,0,1,0,1,1,0,1,1,0,0,1,1},
	{1,1,0,1,1,0,1,1,1,0,1,0,1,1,0,0,0,1,1,0,0,1,0,1,1,1,1,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,1,0,0,1,1,1,0,0,1,1,0,1,0,0,0,0,1,1,1,1},
	{1,1,0,1,0,0,0,1,1,1,1,0,0,1,1,0,1,0,1,1,0,1,1,0,0,0,1,0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,1,1,0,0,1,0,1,0,0,1,0,0,1,1,1,0,1,1,1,1}
};

u64 k[ROUNDS] = { 0 };

std::mt19937 rng;
std::uniform_int_distribution<u64> uni_dist(0x0ull, WORD_MASK);


//////////////////////////////////////////
//          HELPERS
//////////////////////////////////////////
char* binary(u64 x) {
	char *r = new char[BLOCK_SIZE+1];
	for (int i = 0; i < BLOCK_SIZE; ++i)
		r[i] = ((x >> (BLOCK_SIZE-1-i)) & 0x1) ? '1' : '0';
	r[BLOCK_SIZE] = '\0';
	return r;
}

// Hamming weight (64-bit words)
int weight(u64 x) {
	int r = 0;
	while (x) {
		if (x & 0x1)
			r++;
		x >>= 1;
	}
	return r;
}

// longest common subsequence, assumes len(x) = len(y) = BLOCK_SIZE+1 (where the last is \0)
unsigned char lcs(char *x, char *y) {
    unsigned char L[(BLOCK_SIZE+1)][(BLOCK_SIZE+1)];
    for(int i = 0; i <= BLOCK_SIZE; i++) {
        for(int j = 0; j <= BLOCK_SIZE; j++) {
            if(i == 0 || j == 0)
                L[i][j] = 0;
            else if(x[i-1] == y[j-1])
                L[i][j] = L[i-1][j-1] + 1;
            else
                L[i][j] = (L[i-1][j] > L[i][j-1]) ? L[i-1][j] : L[i][j-1];
        }
    }
    return L[BLOCK_SIZE][BLOCK_SIZE];
}

//////////////////////////////////////////
//          CIPHER STUFF
//////////////////////////////////////////
// Rotate. As in the paper, positive amount is left, negative amount is right.
u64 rotate(u64 x, int p) {
	if (p >= BLOCK_SIZE || p <= -BLOCK_SIZE)
		perror("Bad rotation amount!\n");
		
	return (p > 0) ? 
		((x << p) | (x >> (BLOCK_SIZE-p))) & WORD_MASK :
		((x >> (-p)) | (x << (BLOCK_SIZE+p))) & WORD_MASK;

	/* OLD METHOD
	return (p > 0) ? 
		(((x & (((0x1ull << p) - 1) << (BLOCK_SIZE-p))) >> (BLOCK_SIZE-p)) ^ (x << p)) & WORD_MASK : 
		(((x & (((0x1ull << (BLOCK_SIZE+p)) - 1) << (-p))) >> (-p)) ^ (x << (BLOCK_SIZE+p))) & WORD_MASK;
	*/
}

// F function
u64 F(u64 x) {
	return (rotate(x,1) & rotate(x,8)) ^ rotate(x,2);
}

// Toy version of F, rotation by 8 is replaced by p
u64 F_toy(u64 x) {
	return (rotate(x,1) & rotate(x,8)); // ^ rotate(x,2);
}

// Key schedule
void key_schedule() {
	u64 tmp;
	for (int i = KEY_WORDS; i < ROUNDS; ++i) {
		tmp = rotate(k[i-1], -3);
		if (KEY_WORDS == 4)
			tmp ^= k[i-3];
		tmp ^= rotate(tmp, -1);
		k[i] = k[i-KEY_WORDS] ^ z[CONST_J][(i-KEY_WORDS) % 62] ^ tmp ^ CONST_C;
		//k[i] = k[i-KEY_WORDS] ^ tmp;
	}
	if (PRINT_ROUND_KEYS) {
		for (int i = 0; i < ROUNDS; ++i)
			printf("%s\n", binary(k[i]));
		printf("\n\n");
	}
			//printf("k[%2d] : %s   wt : %d\n", i, binary(k[i]), weight(k[i]));
}

// Encryption for num_rounds rounds
void encrypt(u64 &x, u64 &y, int num_rounds) {
	u64 tmp;
	for (int i = 0; i < num_rounds; ++i) {
		tmp = x;
		x = y ^ F(x) ^ k[i];
		y = tmp;
	}
}

// Encrypt for the full number of rounds
void encrypt(u64 &x, u64 &y) {
	encrypt(x, y, ROUNDS);
}

// Decryption for num_rounds rounds
void decrypt(u64 &x, u64 &y, int num_rounds) {
	u64 tmp;
	for (int i = 0; i < num_rounds; ++i) {
		tmp = y;
		y = x ^ F(y) ^ k[ROUNDS-i-1];
		x = tmp;
	}
}

// Decrypt for the full number of rounds
void decrypt(u64 &x, u64 &y) {
	decrypt(x, y, ROUNDS);
}

void setup_random_key() {
	// key schedule
	for (int i = 0; i < KEY_WORDS; ++i)
		k[i] = uni_dist(rng) & WORD_MASK;
	key_schedule();
}

//////////////////////////////////////////
//          CIPHER TESTING
//////////////////////////////////////////
void run_test_vectors() {
	u64 x, y, ex, ey;
	if (BLOCK_SIZE == 16 && KEY_WORDS == 4) {
		k[3] = 0x1918;	k[2] = 0x1110;	k[1] = 0x0908;	k[0] = 0x0100; x = 0x6565;	y = 0x6877;	ex = 0xc69b; ey = 0xe9bb;
	}
	if (BLOCK_SIZE == 24 && KEY_WORDS == 3) {
		k[2] = 0x121110; k[1] = 0x0a0908; k[0] = 0x020100; x = 0x612067; y = 0x6e696c; ex = 0xdae5ac; ey = 0x292cac;
	}
	if (BLOCK_SIZE == 24 && KEY_WORDS == 4) {
		k[3] = 0x1a1918; k[2] = 0x121110; k[1] = 0x0a0908; k[0] = 0x020100; x = 0x726963; y = 0x20646e;	ex = 0x6e06a5; ey = 0xacf156;
	}
	if (BLOCK_SIZE == 32 && KEY_WORDS == 3) {
		k[2] = 0x13121110; k[1] = 0x0b0a0908; k[0] = 0x03020100; x = 0x6f722067; y = 0x6e696c63; ex = 0x5ca2e27f; ey = 0x111a8fc8;		
	}
	if (BLOCK_SIZE == 32 && KEY_WORDS == 4) {
		k[3] = 0x1b1a1918; k[2] = 0x13121110; k[1] = 0x0b0a0908; k[0] = 0x03020100; x = 0x656b696c;	y = 0x20646e75; ex = 0x44c8fc20; ey = 0xb9dfa07a;		
	}
	if (BLOCK_SIZE == 48 && KEY_WORDS == 2) {
		k[1] = 0x0d0c0b0a0908; k[0] = 0x050403020100; x = 0x2072616c6c69; y = 0x702065687420; ex = 0x602807a462b4; ey = 0x69063d8ff082;
	}
	if (BLOCK_SIZE == 48 && KEY_WORDS == 3) {
		k[2] = 0x151413121110; k[1] = 0x0d0c0b0a0908; k[0] = 0x050403020100; x = 0x746168742074; y = 0x73756420666f; ex = 0xecad1c6c451e; ey = 0x3f59c5db1ae9;
	}	
	if (BLOCK_SIZE == 64 && KEY_WORDS == 2) {
		k[1] = 0x0f0e0d0c0b0a0908; k[0] = 0x0706050403020100; x = 0x6373656420737265; y = 0x6c6c657661727420; ex = 0x49681b1e1e54fe3f; ey = 0x65aa832af84e0bbc;
	}
	if (BLOCK_SIZE == 64 && KEY_WORDS == 3) {	
		k[2] = 0x1716151413121110; k[1] = 0x0f0e0d0c0b0a0908; k[0] = 0x0706050403020100; x = 0x206572656874206e; y = 0x6568772065626972; ex = 0xc4ac61effcdc0d4f; ey = 0x6c9c8d6e2597b85b;
	}
	if (BLOCK_SIZE == 64 && KEY_WORDS == 4) {
		k[3] = 0x1f1e1d1c1b1a1918; k[2] = 0x1716151413121110; k[1] = 0x0f0e0d0c0b0a0908; k[0] = 0x0706050403020100; x = 0x74206e69206d6f6f; y = 0x6d69732061207369; ex = 0x8d2b5579afc8a3a0; ey = 0x3bf72a87efe7b868;
	}
	
	
	key_schedule();
	encrypt(x, y);
	if (x != ex || y != ey)
		printf("Test-std::vector mismatch! %016llx %016llx <=> %016llx %016llx\n", x, y, ex, ey);
}

void test_enc() {	
	u64 x,y,a,b;
	
	for (int j = 0; j < 1000000; ++j) {
		// Draw random plaintext
		a = uni_dist(rng) & WORD_MASK;
		b = uni_dist(rng) & WORD_MASK;
		x = a;
		y = b;
		
		// Draw random master key
		for (int i = 0; i < KEY_WORDS; ++i)
			k[i] = uni_dist(rng) & WORD_MASK;
		key_schedule();
		
		// Encrypt + decrypt
		encrypt(x,y);
		decrypt(x,y);
		
		if (x != a || y != b)
			printf("encrypt/decrypt failed!\n");
	}
}

//////////////////////////////////////////
//          DIFFERENTIAL STUFF
//////////////////////////////////////////
void distribution_test() {
	int RUNS = 1000;
	unsigned int weight_count[ROUNDS][17] = { { 0 } };
	
	for (int i = 0; i < RUNS; ++i) {
		k[3] = uni_dist(rng);	
		k[2] = uni_dist(rng);
		k[1] = uni_dist(rng);
		k[0] = uni_dist(rng);
		key_schedule();
		
		/*u64 x = uni_dist(rng);
		u64 y = uni_dist(rng);
		u64 x2 = x ^ 0x1ull;
		u64 y2 = y;
		
		encrypt(x,y,ROUNDS);
		encrypt(x2,y2,ROUNDS);
		*/
		for (int j = 0; j < ROUNDS; ++j)
			weight_count[j][weight(k[j])]++;
	}
	for (int i = 0; i < 17; ++i) {
		for (int j = 0; j < ROUNDS; ++j)
			printf("%f\t", (double)weight_count[j][i] / RUNS);
		printf("\n");
	}
}

#define table(i,j) (table[(i)*(WORD_MASK+1) + (j)])

void diff_dist_table() {
	if (BLOCK_SIZE == 16) {
		//unsigned int table[WORD_MASK + 1][WORD_MASK + 1] = { { 0 } };
		unsigned short *table = (unsigned short*) calloc((WORD_MASK+1)*(WORD_MASK+1), sizeof(unsigned short));
		
		u64 x, y, _tmp;
		for (x = 0; x <= WORD_MASK; ++x)
			for (y = 0; y <= WORD_MASK; ++y)
				table((x^y), (F(x)^F(y)))++;

		u64 best = 0;
		std::vector<std::pair<u64,u64>> best_pairs;
		/*for (x = 0; x <= WORD_MASK; ++x) {
			for (y = 0; y <= WORD_MASK; ++y) {
				if (table(x,y) > best)
					best = table(x,y);
			}
		}
		printf("best entry : %d\n", best);
		return;
		for (x = 1; x <= WORD_MASK; ++x) {
			if (table(x,x) > best) {
				best = table(x,x);
				printf("new best on diagonal %016llx count %d\n", x, table(x,x));
			}
		}
		return;
		*/
		for (x = 1; x <= WORD_MASK; ++x) {
			for (y = 1; y <= WORD_MASK; ++y) {
				_tmp = table(x,y) * table(y,x);
				if (_tmp == best)
					best_pairs.push_back(std::make_pair(x,y));				
				else if (_tmp > best) {
					best_pairs.clear();
					best = _tmp;
					best_pairs.push_back(std::make_pair(x,y));
				}
			}
		}		
		printf("listing best difference pairs found :\n");
		for (std::vector<std::pair<u64,u64>>::iterator it = best_pairs.begin(); it != best_pairs.end(); ++it)
			printf("(a,b) = (%016llx, %016llx) = (%s, %s)    counts : %d and %d\n", (*it).first, (*it).second, binary((*it).first), binary((*it).second), table((*it).first, (*it).second), table((*it).second, (*it).first));
		
		/*best = 0;
		u64 v;
		best_pairs.clear();
		for (x = 0; x <= WORD_MASK; ++x) {
			v = 0;
			for (y = 0; y <= WORD_MASK; ++y)
				v += table(x,y) * table(y,x);
			
			if (v == best)
				best_pairs.push_back(std::make_pair(x,v));
			
			if (v > best) {
				best_pairs.clear();
				best_pairs.push_back(std::make_pair(x,v));
				best = v;
			}
		}
		printf("listing best input differences a found :\n");
		for (std::vector<std::pair<u64,u64>>::iterator it = best_pairs.begin(); it != best_pairs.end(); ++it)
			printf("a = %016llx = %s    count : %d\n", (*it).first, binary((*it).first), (*it).second);
		*/
	}
	else
		perror("Can't do difference distribution table for n > 16!\n");
}

void test_differential(u64 alphaL, u64 alphaR, u64 betaL, u64 betaR, int rounds_to_test) {
	// set up random master key and run key schedule
	setup_random_key();

	u64 x, y, m1L, m1R, m2L, m2R;
	u64 good = 0, counter = 0;
	x = 0;
	//while (x <= WORD_MASK) {
	//	for (i = 0; i < alphaL; ++i) {
	for (x = 0; x <= WORD_MASK; x+=2) {
			for (y = 0; y <= 0xFFF; ++y) {
				m1L = x;
				m1R = y;

				m2L = x ^ alphaL;
				m2R = y ^ alphaR;
				
				encrypt(m1L, m1R, rounds_to_test);
				encrypt(m2L, m2R, rounds_to_test);
				
				if ((m1L^m2L) == betaL && (m1R^m2R) == betaR)
					good++;
				counter++;
			}
	//		x++;
	//	}
	//	x += alphaL;
	}
	printf("experimental differential prob : %u/%u\n", good, counter);
}

// Do iterated differential attack of the form [a 0] -> [b a] -> [0 b] -> [b 0] -> [a b] -> [0 a] -> [a 0]
// a and b are differences of BLOCK_SIZE bits
void diff_attack(u64 alphaL, u64 alphaR, u64 betaL, u64 betaR, int rounds_attacked) {
	setup_random_key();

	std::map<unsigned int, unsigned int> freq;


	unsigned int *counts = (unsigned int*) calloc((WORD_MASK+1), sizeof(unsigned int)); //[WORD_MASK + 1] = { 0 };

	u64 x, y, m1L, m1R, m2L, m2R, keyguess, i;
	clock_t start_time = clock();

	for (y = 0; y <= 0xFFFF; ++y) {
		x = 0;
		while (x <= WORD_MASK) {
			for (i = 0; i < alphaL; ++i) {
				m1L = x;
				m1R = y;

				m2L = x ^ alphaL;
				m2R = y ^ alphaR;
				
				encrypt(m1L, m1R, rounds_attacked);
				encrypt(m2L, m2R, rounds_attacked);
				
				if ( (F(m1R) ^ F(m2R) ^ m1L ^ m2L) == betaL ) {
					for (keyguess = 0; keyguess <= WORD_MASK; ++keyguess) {
						if ( (F(F(m1R)^m1L^keyguess) ^ 
							  F(F(m2R)^m2L^keyguess) ^ m1R ^ m2R) == betaR )
							counts[keyguess]++;
					}
				}					
				x++;
			}
			x += alphaL;
		}
	}

	for (u64 x = 0; x <= WORD_MASK; ++x)
		freq[counts[x]]++;

	printf("allCounts[correct key] = allCounts[%016llx] = %d\n", k[rounds_attacked-1], counts[k[rounds_attacked-1]]);
	printf("Time spent : %f sec.\n", (double)(clock()-start_time)/CLOCKS_PER_SEC);

	for (std::map<unsigned int, unsigned int>::iterator it = freq.begin(); it != freq.end(); ++it)
		printf("%u\t%u\n", (*it).first, (*it).second);

	free(counts);
}

void test_2r_difference(u64 diff) {
	u64 x, c = 0;
	for (x = 0; x <= WORD_MASK; ++x) {
		if ( (F(F(x)) ^ F(F(x^diff))) == diff )
			c++;
	}
	printf("difference %016llx -> * -> %016llx happened %d times out of %d\n", diff, diff, c, (WORD_MASK+1));
}

// determine the diagonal of the DDT in time 2^n memory 2^n
// see the rump session talk from Eurocrypt 2013
void ddt_diagonal() {
	std::map<u64, std::vector<u64>> M;
	u64 x, diff;
	for (x = 0; x <= WORD_MASK; ++x)
		M[(x ^ F(x))].push_back(x);
	
	u64 best = 0, val = 0;
	int len, p, q;
	std::set<u64> good_diffs;
	std::set<u64> seen_diffs;
	
	for (std::map<u64, std::vector<u64>>::iterator it = M.begin(); it != M.end(); ++it) {
		len = ((*it).second).size();
		if (len > 1) {			
			for (p = 0; p < len; ++p) {
				for (q = p+1; q < len; ++q) {
					val = 0;
					diff = ((*it).second)[p]^((*it).second)[q];
					if (!seen_diffs.count(diff)) {
						seen_diffs.insert(diff);
						for (x = 0; x <= WORD_MASK; ++x) {
							if ((F(x)^F(x^diff)) == diff)
								val++;
						}
						if (val >= best) {
							best = val;
							//good_diffs.insert(diff);
							printf("diff %016llx\tcount %d\n", diff, val);
						}
					}
				}
			}
		}
	}
	printf("differences seen : %d\n", seen_diffs.size());
	/*for (set<u64>::iterator it = good_diffs.begin(); it != good_diffs.end(); ++it)
		printf("difference %04llx\n", *it);
	return;
	printf("%d\n", best);
	best = 0;
	u64 d = 0x0000000000003214^0x0000000000009361;
	for (x = 0; x <= WORD_MASK; ++x)
		if ((F(x) ^ F(x^d)) == d)
			best++;
	printf("c : %d\n", best);*/
}

//////////////////////////////////////////
//          ROTATIONAL STUFF
//////////////////////////////////////////
void rotational_approx() {
	u64 x, count, best = 0, rot, Fx;
	u64 *rot_table = (u64 *)calloc(BLOCK_SIZE, sizeof(u64));

	for (x = 0; x <= WORD_MASK; ++x) {
		Fx = F(x);
		for (rot = 0; rot < BLOCK_SIZE; ++rot) {
			if (rotate(Fx,rot) == x)
				rot_table[rot]++;
		}
	}
	printf("== Number of inputs x s.t. F(x <<< alpha) == F(x) <<< alpha\n");
	for (rot = 0; rot < BLOCK_SIZE; ++rot)
		printf("%d\t%d\n", rot, rot_table[rot]);
	
	return;
	
	printf("== Searching for (alpha,e) maximizing Pr(x = (F(x) <<< alpha) + e) ==\n");
	std::vector<std::pair<u64,u64>> best_pairs;
	for (int alpha = 0; alpha < BLOCK_SIZE; ++alpha) {
		for (u64 e = 0; e <= WORD_MASK; ++e) {
			count = 0;
			for (x = 0; x <= WORD_MASK; ++x) {
				if (x == (rotate(F(x), alpha) ^ e))
					count++;
			}
			if (count == best)
				best_pairs.push_back(std::make_pair(alpha,e));
			if (count > best) {
				best = count;
				best_pairs.clear();
				best_pairs.push_back(std::make_pair(alpha,e));
			}
		}
	}
	printf("The best count for x = (F(x) <<< alpha) + e was : %d\nHere is the list:\n", best);
	for (unsigned int i = 0; i < best_pairs.size(); ++i)
		printf("(alpha, e) = (%d, %016llx)\n", best_pairs[i].first, best_pairs[i].second);
	
}

void differences_to_zero_diff() {
	unsigned short *list = (unsigned short*) calloc((WORD_MASK+1), sizeof(unsigned short));

	for (u64 x = 0; x <= WORD_MASK; ++x) {
		for (u64 y = 0; y <= WORD_MASK; ++y) {
			if ((F(x) ^ F(y)) == 0)
				list[(x^y)]++;
		}
	}
	for (u64 i = 0; i <= WORD_MASK; ++i)
		if (list[i] == 128)
			printf("%s\t%d\n", binary(i), list[i]);

	free(list);
}


//////////////////////////////////////////
//          WEAK KEY STUFF
//////////////////////////////////////////
void weak_keys() {
	k[0] = 0xaaaaaaaaaaaaaaaa & WORD_MASK;
	k[1] = 0xaaaaaaaaaaaaaaaa & WORD_MASK;
	k[2] = 0xaaaaaaaaaaaaaaba & WORD_MASK;
	k[3] = 0xaaaaaaaaaaaaaaba & WORD_MASK;
	//k[2] = 0xaaaaaaaaaaaaaaba & WORD_MASK; // (rotate_right(k[0], 3) ^ 0xaaaaaaaaaaaaaaaa) & WORD_MASK;
	//k[3] = 0xaaaaaaaaaaaaaaaa & WORD_MASK; //(rotate_left(k[1], 3) ^ 0xaaaaaaaaaaaaaaaa) & WORD_MASK;
	key_schedule();
}

void test_rotational() {
	// set up random master key and run key schedule
	for (int i = 0; i < KEY_WORDS; ++i)
		k[i] = uni_dist(rng) & WORD_MASK;
	key_schedule();
	
	unsigned short lol[BLOCK_SIZE] = { 0 };
	u64 x,y;	
	for (u64 i = 0; i <= WORD_MASK; ++i) {
		for (u64 j = 0; j <= 0x4000; ++j) {
			x = i;
			y = j;
			
			encrypt(x,y,2);
			
			for (int m = 0; m < BLOCK_SIZE; ++m)
				if (x == rotate(i,m))
					lol[m]++;
		}
	}
	for (int m = 0; m < BLOCK_SIZE; ++m)
		printf("rotation %d : %d\n", m, lol[m]);
}

void test_key_rotation() {
	unsigned int freq[ROUNDS] = { 0 };
	unsigned int RUNS = 0x8000;
	u64 *kcopy = (u64*) malloc(ROUNDS*sizeof(u64));
	
	for (unsigned int j = 0; j < RUNS; ++j) {
		// set up random master key and run key schedule
		for (int i = 0; i < KEY_WORDS; ++i)
			k[i] = uni_dist(rng) & WORD_MASK;
		key_schedule();
		
		// copy round keys		
		memcpy(kcopy, k, ROUNDS*sizeof(u64));
		
		// rotate master key and do key schedule again
		for (int i = 0; i < KEY_WORDS; ++i)
			k[i] = rotate(k[i], 2);
		key_schedule();
		
		// check longest common subsequence of k[i] and (k'[i] >>> 2)
		for (int i = 0; i < ROUNDS; ++i)
			freq[i] += lcs(binary(kcopy[i]), binary(rotate(k[i], -2)));
			//printf("k[%2d] = %s and (k'[%2d] >>> 2) = %s have LCS length %d\n", i, binary(kcopy[i]), i, binary(rotate(k[i], -2)), );
		
	}
	for (int i = 0; i < ROUNDS; ++i)
		printf("%d\t%f\n", i, (double)freq[i] / RUNS);
	
	free(kcopy);
}

void is_F_balanced() {
	unsigned short *abe = (unsigned short*) calloc((WORD_MASK+1), sizeof(unsigned short));
	for (u64 x = 0; x <= WORD_MASK; ++x)
		abe[F(x)]++;
	
	unsigned short t = abe[0];
	for (u64 x = 0; x <= WORD_MASK; ++x) {
		if (abe[x] != 0 && abe[x] != t)
			printf("not balanced!\n");
	}
}

void generate_key_relations(unsigned int rotation) {
	int ROUNDS_TO_DO = 5;

	printf("B.<");
	for (int r = 0; r < ROUNDS_TO_DO; ++r) {
		for (int i = 0; i < BLOCK_SIZE; ++i)
			printf("k%02d%02d,kp%02d%02d,", r,i,r,i);
	}
	printf("> = BooleanPolynomialRing()\n"); 

	// rotational relation between k and kp
	printf("L = [");
	int eidx = 0;
	for (int r = 0; r < ROUNDS_TO_DO; ++r) {
		for (int i = 0; i < BLOCK_SIZE; ++i) {
			printf("kp%02d%02d + k%02d%02d + 1,", r, i, r, (i-rotation)%BLOCK_SIZE);
			eidx++;
		}
	}
	
	// round keys to master key relations
	eidx = 0;
	for (int r = KEY_WORDS; r < ROUNDS_TO_DO; ++r) {
		for (int i = 0; i < BLOCK_SIZE; ++i) {
			switch (i) {
				case 0:
					printf("k%02d%02d + k%02d%02d + k%02d%02d + k%02d%02d", r, i, (r-1), ((i+4)%BLOCK_SIZE), (r-1), ((i+3)%BLOCK_SIZE), (r-KEY_WORDS), i);
					if (!z[CONST_J][(r-KEY_WORDS) % 62])
						printf(" + 1");
					break;
				case 1:
					printf("k%02d%02d + k%02d%02d + k%02d%02d + k%02d%02d + 1", r, i, (r-1), ((i+4)%BLOCK_SIZE), (r-1), ((i+3)%BLOCK_SIZE), (r-KEY_WORDS), i);
					break;
				default:
					printf("k%02d%02d + k%02d%02d + k%02d%02d + k%02d%02d", r, i, (r-1), ((i+4)%BLOCK_SIZE), (r-1), ((i+3)%BLOCK_SIZE), (r-KEY_WORDS), i);
					break;
			}
			if (KEY_WORDS == 4)
				printf(" + k%02d%02d", (r-3), ((i+1)%BLOCK_SIZE));
			printf(",");
			eidx++;
		}
	}
	// round keys to master key relations
	for (int r = KEY_WORDS; r < ROUNDS_TO_DO; ++r) {
		for (int i = 0; i < BLOCK_SIZE; ++i) {
			switch (i) {
				case 0:
					printf("kp%02d%02d + kp%02d%02d + kp%02d%02d + kp%02d%02d", r, i, (r-1), ((i+4)%BLOCK_SIZE), (r-1), ((i+3)%BLOCK_SIZE), (r-KEY_WORDS), i);
					if (!z[CONST_J][(r-KEY_WORDS) % 62])
						printf(" + 1");					
					break;
				case 1:
					printf("kp%02d%02d + kp%02d%02d + kp%02d%02d + kp%02d%02d + 1", r, i, (r-1), ((i+4)%BLOCK_SIZE), (r-1), ((i+3)%BLOCK_SIZE), (r-KEY_WORDS), i);
					break;
				default:
					printf("kp%02d%02d + kp%02d%02d + kp%02d%02d + kp%02d%02d", r, i, (r-1), ((i+4)%BLOCK_SIZE), (r-1), ((i+3)%BLOCK_SIZE), (r-KEY_WORDS), i);
					break;
			}
			if (KEY_WORDS == 4)
				printf(" + kp%02d%02d", (r-3), ((i+1)%BLOCK_SIZE));
			printf(",");
			eidx++;
		}
	}
	printf("]\n");
}

void key_difference() {
	u64 *kcopy = (u64*) malloc(ROUNDS*sizeof(u64));
	for (int i = 0; i < KEY_WORDS; ++i)
		k[i] = uni_dist(rng) & WORD_MASK;
	key_schedule();
	memcpy(kcopy, k, ROUNDS*sizeof(u64));
	
	k[0] = kcopy[0];
	k[1] = kcopy[1]^0x8000000000000000;
	key_schedule();
	
	for (int i = 0; i < ROUNDS; ++i)
		printf("%s\n", binary(k[i]^kcopy[i]));
}

//////////////////////////////////////////
//     IMPOSSIBLE DIFFERENTIAL STUFF
//////////////////////////////////////////
std::string rotchar(std::string x, unsigned int p) {
	std::string ret(BLOCK_SIZE, '0');
	for (int i = 0; i < BLOCK_SIZE; ++i)
		ret[i] = x[(i+p) % BLOCK_SIZE];
	return ret;
}

unsigned short cweight(std::string x) {
	unsigned short c = 0;
	for (int i = 0; i < BLOCK_SIZE; ++i)
		c = x[i] == '0' ? c+1 : c;
	return BLOCK_SIZE-c;
}

std::string cxor(std::string a, std::string b) {
	std::string x(BLOCK_SIZE, '0');
	for (int i = 0; i < BLOCK_SIZE; ++i) {
		if (a[i] == '*' || b[i] == '*')
			x[i] = '*';
		else
			x[i] = (char)(48 + ((a[i] - 48) ^ (b[i] - 48)));
	}
	return x;
}

std::string rot_func(std::string in) {
	std::string ret(BLOCK_SIZE, '0');
	
	for (int i = 0; i < BLOCK_SIZE; ++i) {
		if (in[i] == '*') {
			ret[(i-1+BLOCK_SIZE)%BLOCK_SIZE] = '*';
			ret[(i-2+BLOCK_SIZE)%BLOCK_SIZE] = '*';
			ret[(i-8+BLOCK_SIZE)%BLOCK_SIZE] = '*';
		}
		else if (in[i] == '1') {
			ret[(i-1+BLOCK_SIZE)%BLOCK_SIZE] = '*';
			ret[(i-8+BLOCK_SIZE)%BLOCK_SIZE] = '*';
				
			if (ret[(i-2+BLOCK_SIZE)%BLOCK_SIZE] == '0')
				ret[(i-2+BLOCK_SIZE)%BLOCK_SIZE] = '1';
		}
	}
	return ret;
}

void impossible_diff_attack() {
	// out pattern: **************** 0******1******0*
	// in pattern : **************** 0******1******0*

	setup_random_key();

	// make set with all key candidates
	std::set<u64> key_candidates;
	for (u64 kg = 0; kg <= WORD_MASK; ++kg)
		key_candidates.insert(kg);
	
	// start timing
	//omp_set_num_threads(THREADS);
	clock_t start = clock();
	double start_d = omp_get_wtime();

	/*#pragma omp parallel default(none), shared(key_candidates), private(k)
	{
		#pragma omp single
		{
			printf("# of threads : %d\n", omp_get_num_threads());
		}
		int tid = omp_get_thread_num();
		printf("thread %d reporting for duty!\n", tid);
	*/
		u64 x, y, aL, aR, bL, bR, aQ, bQ, Fx, Fxp;
		u64 ind, outd1, outd2;
		std::set<u64>::iterator candit;
		unsigned int removedByMe = 0;

		ind = 0x1;
		outd1 = rotate(ind, 7);
		outd2 = rotate(ind, 9);
		
	//	#pragma omp for
		for (x = 0; x <= WORD_MASK; ++x) {
			Fx = F(x);
			
		/*	for (int rot = 0; rot < BLOCK_SIZE; ++rot) {
				ind 	= rotate(0x1, rot);
				outd1 	= rotate(ind, 7);
				outd2 	= rotate(ind, 9);
		*/		
				Fxp = F(x ^ ind);
			
				for (y = 0; y < 0x4000; ++y) {
					// construct input
					aL = x;
					bL = x ^ ind;
					aR = y;
					bR = y ^ Fx ^ Fxp;
					
					// encrypt
					encrypt(aL, aR, 14);
					encrypt(bL, bR, 14);
					
					/*if ( ((aR^bR) & rotate(0x8102, rot)) != rotate(0x0100, rot) )
						printf("shit hit the fan!\n");*/
					
					// key recovery
					aQ = F(aR);
					bQ = F(bR);
					if ( (aQ ^ bQ ^ aL ^ bL) == outd1 || (aQ ^ bQ ^ aL ^ bL) == outd2 ) { // first filter
						//for (candit = key_candidates.begin(); candit != key_candidates.end(); ++candit) {
						//if ( (F(aL^aQ^(*candit)) ^ F(bL^bQ^(*candit)) ^ aR ^ bR) == 0x0 ) {
						for (u64 keyguess = 0; keyguess <= WORD_MASK; ++keyguess) {							
							if ( (F(aL^aQ^keyguess) ^ F(bL^bQ^keyguess) ^ aR ^ bR) == 0x0 ) {
								key_candidates.erase(keyguess);
								removedByMe++;
								//printf("%d ", key_candidates.size()); //    contained : %d\n", key_candidates.size(), key_candidates.count(k[12]));
							}
						}
					}
				}
			//}
		}

		printf("I removed %d key candidates\n", removedByMe);
		//printf("thread %d removed %d keys\n", tid, removedByMe);

	//} // end omp parallel

	/*for (int tid = 0; tid < THREADS; ++tid) {
		for (set<u64>::iterator keyguess = removed_keys[tid].begin(); keyguess != removed_keys[tid].end(); ++keyguess)
			key_candidates.erase((*keyguess));
	}*/

	printf("time.h : %f sec.\nomp.h  : %f sec.\n", (double)(clock()-start)/CLOCKS_PER_SEC, (omp_get_wtime()-start_d));	
	printf("real sub key : %016llx is contained ? %d\nthere are %d possible keys left\n", k[13], key_candidates.count(k[13]), key_candidates.size());
}

std::string dL, dR, tm;
std::string zeroes(BLOCK_SIZE, '0');

int impossible_diff(std::string leftDiff, std::string rightDiff) {
	dL = leftDiff; //(BLOCK_SIZE, '0'); dL[BLOCK_SIZE-1] = '1';
	dR = rightDiff;
	tm = zeroes;
	
	short c = 0;
	do {
		printf("%s %s\n", dL.c_str(), dR.c_str());
		tm = dL;
		dL = cxor(rot_func(dL), dR);
		dR = tm;
		c++;
	} while(cweight(dL) < BLOCK_SIZE || cweight(dR) < BLOCK_SIZE);
	printf("rounds = %d\n", c-1);
	return c-1;
}


void imp_diff_attack3() {
	setup_random_key();

	std::set<u64> remaining_keys;
	for (u64 keyguess = 0; keyguess <= WORD_MASK; ++keyguess)
		remaining_keys.insert(keyguess);

	clock_t start_time = clock();

	std::set<u64>::iterator siter;
	u64 keyguess, x, y, xPlusAlpha, afterFdiff, tmp1, tmp2, m1L, m1R, m2L, m2R, alpha, outdiff1, outdiff2;
	unsigned short rot, xinnerl;
	
	for (rot = 0; rot < BLOCK_SIZE; ++rot) {
		alpha = rotate(0x1, rot);
		outdiff1 = rotate(alpha, 7);
		outdiff2 = rotate(alpha, 9);

		x = 0;
		while (x <= 18860) {
			// use the next 2^rot values of x
			for (xinnerl = 0; xinnerl < alpha; ++xinnerl) {
				afterFdiff = F(x) ^ F(x ^ alpha);
				xPlusAlpha = x ^ alpha;

				for (y = 0; y <= WORD_MASK; ++y) {
					m1L = x;
					m1R = y;

					m2L = xPlusAlpha;
					m2R = y ^ afterFdiff;

					encrypt(m1L, m1R, 14);
					encrypt(m2L, m2R, 14);

					// SNIPPET FOR TESTING THE IMPOSSIBLE 10-ROUND DIFFERENTIAL
					// UNCOMMENT TO USE, AND MAKE SURE USING RIGHT INPUT DIFFERENCE AND ENCRYPTING 10 ROUNDS ABOVE
					/*if ((m1R ^ m2R) == 0x0) {
						if ((m1L ^ m2L) == outdiff1 || (m1L ^ m2L) == outdiff2)
							printf("the impossible differential was seen :(!\n");
					}*/

					if ( (F(m1R) ^ F(m2R) ^ m1L ^ m2L) == outdiff1 || (F(m1R) ^ F(m2R) ^ m1L ^ m2L) == outdiff2) {
						for (siter = remaining_keys.begin(); siter != remaining_keys.end(); ++siter) {
							keyguess = *siter;

							tmp1 = F(m1R) ^ m1L ^ keyguess;
							tmp2 = F(m2R) ^ m2L ^ keyguess;

							if ( (F(tmp1) ^ F(tmp2) ^ m1R ^ m2R) == 0x0 )
								remaining_keys.erase(keyguess);
						}
					}
				}
				x++;
			}
			// skip the next 2^rot values of x
			x += alpha;
		}
		printf("Remaining after rot = %d   :   %d\n", rot, remaining_keys.size());
	}

	printf("Attacked key       : %04x\n", k[13]);
	printf("Remaining keys     : %d    correct key contained (y/n) : %d\n", remaining_keys.size(), remaining_keys.count(k[13]));
	printf("Time (one thread)  : %f s.\n", ((double)clock() - start_time) / CLOCKS_PER_SEC);
}


//////////////////////////////////////////
//     DFS DIFFERENTIAL SEARCH STUFF
//////////////////////////////////////////

void H(u64 diff, u64 &m, u64 &w) {
	u64 m1 = 0;
	u64 m8 = 0;
	w = 0;
	int i;

	for (i = 0; i < BLOCK_SIZE; ++i) {
		if ((diff >> i) & 1) {
			// saw a 1 on position i
			// put zeroes on positions i+1 and i+8
			m1 |= (0x1 << ((i+1)%BLOCK_SIZE));
			m8 |= (0x1 << ((i+8)%BLOCK_SIZE));
		}	
	}
	m = (m1 | m8) ^ WORD_MASK;

	for (i = 0; i < BLOCK_SIZE; ++i) {
		if (
			((m1 >> i) & 1) == 0 &&
			((m8 >> i) & 1) == 1 &&
			((m1 >> ((i+7) % BLOCK_SIZE)) & 1) == 1 &&
			((m8 >> ((i+7) % BLOCK_SIZE)) & 1) == 0
		)
			w ^= (0x1ull << ((i+7) % BLOCK_SIZE));			
	}

	m |= w;
}

const int SEARCH_ROUNDS = 12;
std::map<std::pair<u64,u64>, double> out_diffs;

std::pair<u64,u64> best_outdiff;
double best_prob = 0.0f;
double round_best_prob[SEARCH_ROUNDS+1] = { 0.0f };

std::map<std::pair<u64,u64>, std::map<int, u64> > prob_counts;

void diff_BB(u64 leftDiff, u64 rightDiff, int round, double prob) {
	round_best_prob[round] = (prob > round_best_prob[round]) ? prob : round_best_prob[round];

	if (round == SEARCH_ROUNDS) {
		// we reached a leaf, what is the probability??
		std::pair<u64,u64> diff = std::make_pair(leftDiff, rightDiff);
		
		//if (leftDiff == 0x0100 && rightDiff == 0x0000)
		prob_counts[diff][floor(log2(prob))]++;

 		out_diffs[diff] += prob;
 		if (out_diffs[diff] > best_prob) {
 			best_prob = out_diffs[diff];
 			printf("New best %016llx %016llx   prob. : 2^{%f}\n", diff.first, diff.second, log2(best_prob));
 			best_outdiff = diff;
 		}
		return;
	}

	u64 m, w, x, beta, rotLeftDiff = rotate(leftDiff, 2);
	const u64 v = 0;
	double next_prob;
	unsigned int c, branches;

	H(leftDiff, m, w);
	branches = (1 << (BLOCK_SIZE - weight(m)));
	next_prob = prob * (1.0 / (double)branches);

	//if (branches > 128) // branches is always a power of 2
	//	return;

	if (next_prob >= (round_best_prob[round+1] / 4096.0f) ) {
		x = 0;
		c = 0;
		while (c < branches && log2(best_prob) < -29.5) {
		    if ((x & m) == v) {
		    	beta = (x | (rotate(x, 7) & w)) ^ rotLeftDiff;
				
				diff_BB((beta^rightDiff), leftDiff, round+1, next_prob);
	    		//branch_best = tmp_p < branch_best ? tmp_p : branch_best;
	    	
		        x++;
		        c++;
		    }
		    else
		        x += (x & m) ^ v;
		}
	}
	return;
}

void test_H() {
	std::set<u64> actual_diffs, method_diffs;
	int ok_flag = 1, out_betas, c;
	u64 x, alpha, beta, m, v, w;

	while (ok_flag) {
		alpha = uni_dist(rng); //0x0007;
		actual_diffs.clear();
		method_diffs.clear();

		// real
		for (x = 0; x <= WORD_MASK; ++x) {
			beta = F(x) ^ F(x ^ alpha);
			actual_diffs.insert(beta);
		}

		H(alpha, m, w);
		v = 0;

		out_betas = (1 << (BLOCK_SIZE - weight(m)));
		printf("alpha    : %s\n", binary(alpha));
		printf("m        : %s\n", binary(m));
		printf("w        : %s\n", binary(w));
		printf("v        : %s\n", binary(v));
		printf("branches : %d\n", out_betas);

		// alternative method
		x = 0;
		c = 0;
		while (c < out_betas) {			
		    if ((x & m) == v) {
		    	beta = (x | (rotate(x, 7) & w)) ^ rotate(alpha, 2);
		    	method_diffs.insert(beta);
		    	
		        x++;
		        c++;
		    }
		    else
		        x += (x & m) ^ v;
		}
		printf("done\n");


		// first compare sizes
		if (actual_diffs.size() != method_diffs.size()) {
			printf("not same sizes\n");
			ok_flag = 0;
		}
		else {
			// if sizes ok, compare elementwise
			for (std::set<u64>::iterator it = actual_diffs.begin(); it != actual_diffs.end(); ++it) {
				if (method_diffs.count(*it) == 0) {
					printf("element not contained in ggen\n");
					ok_flag = 0;
					break;
				}
			}
		}

		if (!ok_flag) {
			printf("i_diff = %s\n-------\n", binary(alpha));
			std::set<u64>::iterator jt = method_diffs.begin();
			for (std::set<u64>::iterator it = actual_diffs.begin(); it != actual_diffs.end(); ++it) {
				printf("%s\n", binary(*it ^ *jt));
				jt++;
			}
			return;
		}
	}
}

void check_diff_effect() {
	std::pair<u64,u64> beta;
	for (std::map<std::pair<u64,u64>, std::map<int, u64>>::iterator it = prob_counts.begin(); it != prob_counts.end(); ++it) {
		beta = it->first;

		if (log2(out_diffs[beta]) >= -33) {
				printf("beta = %016llx %016llx\n", beta.first, beta.second);

			for (std::map<int, u64>::iterator jt = prob_counts[beta].begin(); jt != prob_counts[beta].end(); ++jt)
				printf("%d\t%u\n", jt->first, jt->second);

			printf("total p = %f\n", log2(out_diffs[beta]));
			printf("\n");
		}
	}
}

int main() {
	rng.seed(time(NULL)); // seed marsenne twister rng

	diff_BB(0x0001, 0x0000, 0, 1.0);
	
	//test_differential(0x0001, 0x0000, 0x1501, 0x0404, 6);
	//diff_attack        (0x0001, 0x0000, 0x0201, 0x0100, 13);
	
	//std::string a(BLOCK_SIZE, '0'); //a[BLOCK_SIZE-1] = '1';
	//std::string b(BLOCK_SIZE, '0'); b[BLOCK_SIZE-8] = '1';
	//impossible_diff(a, b);
	//impossible_diff_attack();
	
	//imp_diff_attack3();

	//key_difference();
	
	//generate_key_relations(2);

	//differences_to_zero_diff();

	//test_differential(0x5555);
	//iterated_diff_attack(0x1111);
	//diff_dist_table();
	//weak_keys();
	//ddt_diagonal();
	//is_F_balanced();
	
	//test_key_rotation();
	//rotational_approx();
	//test_rotational();
	//run_test_vectors();
	//test_enc();

	return 0;
}