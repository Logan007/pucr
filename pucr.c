// #include "pucr.h"
// #include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>

#include <limits.h> // rotl

#include <time.h>

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

static clock_t start, end;
static double cpu_time_used;

#define NO_MAX_GAMMA 0xFF

#define MAX_RLE_RANKS 31

#define LZ	0x0003
#define LIT	0x0002
#define RLE	0x0004
#define MARKED	0x0080

#define LZ_LEN_MAX_GAMMA 15
#define MAX_ESC 9
#define MAX_MTF 9

#define FAST_LANE 0

#define OUT_SIZE 65536

#define DECODE_ITERATIONS 10000

// pthread_mutex_t mutex_mtf;

struct pucr_node {
	uint16_t	rle_count;

	uint16_t	lz_count;
	uint16_t	lz_max_offset;
	uint16_t	lz_cur_offset;
	uint32_t	lz_cost;

	uint16_t	dict_min_offset;
	uint16_t	dict_count;

	uint16_t	seen_before;
	uint16_t	see_next;

	uint16_t	way_to_go;
	uint16_t	next_node;
	uint32_t	bits_to_end;
	uint8_t 	new_esc;
	uint8_t		cur_lit;
};


// --== GENERAL ==--

// DONE: implement speed checks using CPU ticks
// TODO: streamline datatypes
// TODO: handle return value of 'init' function

// --== SPEED ==--

// DONE: allow selection of a simpler alternative which uses sane
//	 in_len-dependant defaults  gaining speed by not searching
//	 for the optimum settings
// TODO: check opportunities to speed up string matching,
//	 keep in mind especially short block sizes:
//	 O(n³) might be less less expensive than O(n) for small n,
//	 do the original hash tables help?
// TODO: make data structures thread safe, especially outbuffer related
// TODO: implement pthread support for linux & co.
//	 definetly for optimize_path which requires more memory then,
//	 maybe for string & rle matching
// DONE: rework int_log2 -- maybe include LUT for 0x0000 <= x <= 0xFFFF
// TODO: squeeze formulas after prooved correctness
// TODO: RLE ranking seems to be expensive performancewise, its effectiveness
//	 should be inspected

// --== FUNCTION ==--

// TODO: check if there is a need for max E. Gamma length mechanism, maybe
//	 implemented seperately for lz count, lz offset, rle rank pointer,
//	 rle length; put according information in header
// TODO: implement an optimization step to find out if rle rank table makes sense
//	 for that particular block to compress, i.e. if more bites are saved if
//	 not using a rle rank table (if no_esc > 1 and '1111' is indicating unranked
//	 rle character, common rle runs of 2 would take 16 or more bits and thus
//	 unfortunately too expensive)
// DONE: implement a shorter 'code' for LZ2, i.e. 'ESC00' instead of 'ESC000'.
//	 at least as long as 'ESC001' is not used, would save one bit per LZ2

void start_clock(void) {
	start = clock();
}


double print_clock(void) {
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	fprintf (stderr,"CPU TIME USED: %f\n", cpu_time_used);
	return (cpu_time_used);
}


// https://hbfs.wordpress.com/2008/08/05/branchless-equivalents-of-simple-functions/
// Steven Pigeon
inline unsigned int sign_extension (int x)
 {
  return x >> (CHAR_BIT*sizeof(int)-1);
 }

/* inline int sign_extension (int32_t x) // only works on LE machine !!!
 {
  union
   {
    // let us suppose long is twice as wide as int
    int64_t w;

    // should be hi,lo on a big endian machine
    struct { int32_t lo, hi; }
   } z = { .w=x };

  return z.hi; // lo on BE machine !!!
 }*/

/* int32_t sign_extension (int32_t x)
 {
  return ~((x<0)-1);
 }*/


inline int32_t abs_bl(int32_t x)
{
  return (x + sign_extension(x)) ^ sign_extension(x);
}

inline int32_t sign(int32_t x) {
    return sign_extension(x) - sign_extension(-x);
}

inline int min(int a, int b)
 {
  return b + ((a-b) & sign_extension(a-b));
 }
inline int max(int a, int b)
 {
  return a + ((b-a) & ~sign_extension(b-a));
 }


static uint8_t int_log2_lut [65536];
static int8_t lz_offset_encoded_length[16][16];

uint32_t pucr_init (void) {

// pthread_mutex_init(&mutex_mtf, NULL);

	// initialize LUT for log2
	uint16_t i;
	int16_t v = -1;
	uint16_t mask = 0xFFFF;

	int_log2_lut[0] = 0; // special definition here
	for (i=1; i < 65535; i ++) {
		if ( (i & mask) != 0 ) {
			mask <<= 1;
			v++;
		}
	int_log2_lut[i] = v;
	}

	// precalculate LUT for lz-bitcount (p = no of real bits) for different lsb-msb length (m = no of pure LSB bits) encodings
	for (int m=0;m<16;m++)  {
		for (int p=0; p<16; p++) {
			lz_offset_encoded_length[m][p]=m + ((p<m?0:p-m+1))*2+1 ; // !!! easier to understand, could be shortened for performance
		}
	}
}

// #define int_log2(x) (x?31-__builtin_clz(x):0)
static uint8_t inline int_log2 (uint16_t x) { // uint16_t x

//	if (x < 65535) {
		return (int_log2_lut[x]);
//	} else {
		/* taken from Henry S. Warren Jr.'s 2nd edition of Hacker's Delight */
/*		uint8_t n = 0;
		if (x <= 0x0000FFFF) { n = n + 16; x = x << 16;}
		if (x <= 0x00FFFFFF) { n = n +  8; x = x <<  8;}
		if (x <= 0x0FFFFFFF) { n = n +  4; x = x <<  4;}
		if (x <= 0x3FFFFFFF) { n = n +  2; x = x <<  2;}
		n = n - (x >> 31);
		return (30 - n);
	}*/
}


// ---=== MOVE TO FRONT   MOVE TO FRONT   MOVE TO FRONT   MOVE TO FRONT ===---


static inline uint8_t move_to_front_encode (uint8_t in_char,
                                     uint8_t alphabet[],
				     uint8_t second_line) {
        // search letter in alphabet
        uint16_t j= 0;
        for (;; j++) if (unlikely(alphabet[j]==in_char)) break;
	uint8_t new = (j >> second_line);
        // move letter more towards the front of the alphabet
        for (uint8_t k=j; k>new; k--) alphabet[k]=alphabet[k-1];
        alphabet[new]=in_char;
        return(j);
}


static inline uint8_t move_to_front_decode (uint8_t in_char,
                                     uint8_t alphabet[],
				     uint8_t second_line) {

        uint8_t ret = alphabet[in_char];
// pthread_mutex_lock (&mutex_mtf);
        // move letter more towards the front of the alphabet
	uint8_t new = (in_char >> second_line);
        for (uint8_t k=in_char; k>new; k--) alphabet[k]=alphabet[k-1];
        alphabet[new]=ret;
// pthread_mutex_unlock (&mutex_mtf);
        return(ret);
}


// ---=== LIT LIT LIT LIT LIT LIT LIT LIT LIT LIT LIT LIT LIT LIT LIT ===---


static uint32_t update_lit_occ_from_graph (struct pucr_node graph[], uint16_t in_len,
		   uint16_t lit_occ[]) {

	uint16_t i;
	for (i=0; i<256; i++) lit_occ[i] = 0;

	for (i=0; i < in_len; i=graph[i].next_node)
		if (graph[i].way_to_go == LIT) {
			lit_occ[graph[i].cur_lit]++;
		}
	return (0);
}


// does not give exact result as it does not take the 'divides' into account, TODO !!!
static int32_t huffiness (uint16_t occ[], uint16_t first, uint16_t last, uint8_t *opt_level, uint8_t cost, uint8_t lit_cost[]) {
	// how 'huffy' are the values? does it make sense to encode them with an implicit huffman tree?

	uint16_t i;
	int32_t sum = 0;

	if (cost) for (i=0; i < 256; i++) lit_cost[i] = 8;

	uint8_t huff_level = 0;
	int32_t max_sum = sum;
	uint16_t q;
	uint16_t j;

	for (i=last-first+1; i>=4; i=q) {
		q = i >> 2;
		for (j=first;j<(first+q);j++) { sum += occ[j]; if (cost) lit_cost[j]--; }// one bit less
		for (j=first+(q<<1);j<first+i;j++) { sum -= occ[j]; if (cost) lit_cost[j]++; } // one bit more
		if (sum > max_sum) {
			max_sum = sum;
			huff_level++;
		} else {
			break; // no need to go any further
		}
	}

	*opt_level = huff_level;
	return (max_sum);
}


// ---=== RLE RLE RLE RLE RLE RLE RLE RLE RLE RLE RLE RLE RLE RLE RLE RLE ===---

static uint32_t generate_rle_ranks (uint16_t rle_occ[], uint8_t rle_rank[], uint8_t *rle_rank_len) {

	uint16_t max_index, max_value;
	int16_t i;

	for (i=0;i<MAX_RLE_RANKS;i++) {
		/* find characters with highest number of occurences
		 * but not below 2 for ranks 3 and following */
		max_value =  (i < 3)?1:(i < 16)?2:4;
		max_index = 0xFFFF;
		for (uint16_t j=0; j < 256 ; j++) {
			if (rle_occ[j] > max_value) {
				/* check that candidate J is not already in list */
				uint16_t k;
				for (k=0;k<i;k++) if ( j==rle_rank[k] ) break;
				if (k == i) {
					max_value = rle_occ[j];
					max_index = j;
				}
			}
		}
		if ( max_index != 0xFFFF ) {
			rle_rank[i] = max_index;
			*rle_rank_len = *rle_rank_len + 1;
		} else {
			break;
		}
	}


// find optimal length for ranked RLE table
uint32_t min = 0; // 0xFFFFFFFF;
for (uint16_t g=0;g<i;g++) min += rle_occ[rle_rank[g]] * 8;
uint16_t idx = 0xFFFF;
for (uint16_t ranked=0;ranked<i;ranked++) {
uint32_t sum = 0;
for (uint16_t g=0;g<ranked;g++) sum += 8 + rle_occ[rle_rank[g]] * (int_log2(g+1)*2+1);
for (uint16_t g=ranked;g<i;g++) sum += rle_occ[rle_rank[g]] * ((ranked==0?0:int_log2(ranked)+1) +8);
if (sum < min) {
min=sum;
idx=ranked;
}
}
 *rle_rank_len = idx+1;

//	*rle_rank_len = i;
	return (0);
}


static uint32_t update_rle_occ_from_graph (struct pucr_node graph[], const uint8_t *inbuf, uint16_t in_len,
		   uint16_t rle_occ[]) {

	uint16_t i;
	for (i=0; i<256; i++) rle_occ[i] = 0;

	for (i=0; i < in_len; i=graph[i].next_node)
		if (graph[i].way_to_go == RLE)
			rle_occ[inbuf[i]]++;

	return (0);
}


static uint32_t count_rle_and_generate_occ_from_inbuf (struct pucr_node graph[], const uint8_t *inbuf, uint16_t in_len,
			   uint16_t rle_occ[]) {

	uint8_t character = inbuf [in_len-1];
	uint16_t rle_count = 0;
	int32_t i; // required (as file size is 16 bit and i may fall below 0); also needed outside loop -- really !!! ???

	for (i=255; i>=0; i--) rle_occ[i] = 0;

	for (i=in_len-1; i >= 0; i--) {

		if (inbuf[i] == character) {
			rle_count++;
		} else {

			// !!! the following line are for creating occurances table, for optimizing later
			/* count number of occuring (real) RLE runs (>1) per character */
			if (rle_count > 1) rle_occ[inbuf[i+1]]++;
			rle_count = 1;
			character = inbuf[i];
		}
		graph[i].rle_count = rle_count;
	}
	// first position
	if (rle_count > 1) rle_occ[inbuf[i+1]]++;
}


// ---=== LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ ===---


static uint32_t find_closer_match (struct pucr_node graph[], const uint8_t *inbuf, uint16_t in_len,
				   uint16_t pattern_start, uint16_t pattern_length, uint16_t min_match_start) {

        uint16_t j,k ;

	// we will be able to skip at least 2 characters when comparing
	// because we know we already have a 2-byte match.
	// even more if rle is longer in that position !!!
	uint16_t off = (graph[pattern_start].rle_count -1);
	off = off<=2?2:off;
// off=0;
//	graph[pattern_start].lz_cur_offset = graph[pattern_start].lz_max_offset; // 0xFFFF;
        for (j=graph[pattern_start].seen_before; (j >= min_match_start) && (j != 0xFFFF); j=graph[j].seen_before) { // check at offset J for a match ...
        	if (graph[pattern_start].rle_count != graph[j].rle_count) continue;
               	k = j + off;
                for (k=j+off; (k-j)<pattern_length; k++) { // ... of length (K-J)
                	if (inbuf[k] != inbuf[pattern_start+(k-j)]) break; // no match
	        }
        	if (k > pattern_start) k = j; // even beyond pattern
	        k -= j; // match length
	        if (k >= pattern_length) {
	                // treatment of a match found
	        	graph[pattern_start].lz_cur_offset = pattern_start-j; // store the 'real' offset
	                break;
                }
        }
}


static uint32_t prepare_last_seen_from_graph (struct pucr_node graph[], const uint8_t *inbuf, uint16_t in_len) {

	uint16_t last_occ [65536];
	memset (last_occ, 0xFF, 0x20000);

	for (uint16_t i=0; i < in_len;i++ ) {
// if (inbuf[i] == inbuf[i+1]) continue; // !!! fast_lane candidate, 
		uint16_t idx = (inbuf[i] << 8) | inbuf[i+1];
		graph[i].seen_before = last_occ[idx];
		last_occ[idx] = i;
 if (inbuf[i] == inbuf[i+1]) continue; // !!! fast_lane candidate, 
	}
}


static uint32_t find_matches (struct pucr_node graph[], const uint8_t *inbuf, uint16_t in_len,
				uint16_t first_pos, uint16_t last_pos,
				uint8_t recursive, uint8_t fast_lane, uint8_t lz_length_max_gamma, uint8_t lz2_bits, uint8_t lz_lsb) {

	// a block should not be much shorter than 256...512 bytes
	if ( (recursive > 0) && ((last_pos - first_pos) > 512) ) {
		uint16_t sep = (last_pos - first_pos) * 43 / 64;
		find_matches (graph, inbuf, in_len, first_pos, first_pos+sep, recursive-1, fast_lane, lz_length_max_gamma, lz2_bits, lz_lsb);
		find_matches (graph, inbuf, in_len, first_pos+sep+1, last_pos, recursive-1, fast_lane, lz_length_max_gamma, lz2_bits, lz_lsb);
	} else {
uint64_t *pattern;
uint64_t *match;

		uint16_t i,j,k, off;
//		for (i=first_pos; i <= last_pos; i +=fast_lane?graph[i].rle_count:1) {// i+=graph[i].rle_count) {
//		for (i=first_pos; i <= last_pos; i+=graph[i].rle_count) {
		for (i=first_pos; i <= last_pos; i+=(fast_lane?graph[i].rle_count:1) ) {

			graph[i].lz_max_offset = 0xFFFF; // !!! neccessary for msb-lsb-offset calculation of number of bits
			graph[i].lz_count =  0x0000; // = no  match yet
graph[i].lz_cost = 0;
			// we will be able to skip at least 2 characters when comparing
			// because we know we already have a 2-byte match.
			// even more if rle is longer in that position
			off = (graph[i].rle_count -1);
			off = off<=2?2:off;

// off = 0;
// uint16_t cnt = 2;

// pattern = &inbuf[i+off];

			for (j=graph[i].seen_before; (uint16_t)(j+1) ; j=graph[j].seen_before) { // check at offset J for a match ...

// cnt++;
				if  (graph[i].rle_count != graph[j].rle_count) continue;
				k = j + off;
				// if a match crosses last_pos, "< in_len" would be advantageous to "<=last_pos" ... !!! may be consider for fast_lane


//				for (k = j+off; (k < i) && (i+(k-j)) <in_len ; k+=4) { // ... of length (K-J)
				for (; (k < i) && ((i+(k-j))<in_len) && ( (k-j) < (3<<(lz_length_max_gamma-1)) ); k++) { // ... of length (K-J)
					if (inbuf[k] != inbuf[i+(k-j)]) break;
// match = &inbuf[k];
// if (*pattern ^ *match) break;
				}

				if (k <= i) {  // we always have a match of at least 2

					k -= j; // match length
//k -= (j+ int_log2(notbe64()>>3))
					// treatment of a match found
					// naive weighting function: the longer the better // replace with weighting/boundary conditions later on !!!
// experimental: try to calculate cost.
// WHY does it result in a compression gain?!


uint32_t gain = (k==2)?16-(3+(lz2_bits>3)?lz2_bits-1:lz2_bits): // '-1' includes some implied huffyness gains
////	               k*8-(1 + int_log2(k-1)*2+1 + 10 + int_log2((i-j-3)>>10)*2+1);
//	               k*8-(1 + int_log2(k-1)*2+1 + lz_lsb + int_log2((i-j-3)>>lz_lsb)*2+1);
//	               k*8-(1 + int_log2(k-1)*2+1 + (i-j-k)<1024?11:10+int_log2((i-j-k)>>10)*2+1);
	               k*8-(1 + int_log2(k-1)*2+1 + (i-j-3)<1024?11:10+int_log2((i-j-3)>>10)*2+1);
//	               k*8-(1 + int_log2(k-1)*2+1 + (i-graph[i].seen_before-3)<1024?11:10+int_log2((i-graph[i].seen_before-3)>>10)*2+1);

// the ranked could be taken care of? the dynamic max-gamma could be taken care of? !!!

//					if (k > graph[i].lz_count) {
					if (gain > graph[i].lz_cost) {
//graph[i].lz_max_offset = i-graph[i].seen_before; // store the 'real' offset
// graph[i].lz_max_offset = cnt; // store the 'real' offset
							graph[i].lz_max_offset = i-j; // store the 'real' offset
							graph[i].lz_count=k; // ; -(off2-i); // k
graph[i].lz_cost = gain;

					}
				}
//				if ( k >= (in_len >> (7+fast_lane))) break; // long enough
			}
		}
fprintf (stderr, "--- FINDING MATCHES %5u ... %5u: ",first_pos,last_pos);
print_clock(); start_clock();
	}
}


struct lz_offset_entry {
	uint16_t offset;
	uint16_t count;
	uint16_t bits_used;
};


int compare_lz_offset( const void * a, const void * b) {

    struct lz_offset_entry *entryA = (struct lz_offset_entry *)a;
    struct lz_offset_entry *entryB = (struct lz_offset_entry *)b;

    return ( (entryB->count*entryB->bits_used) - (entryA->count*entryA->bits_used) );
}


static uint32_t find_optimal_lz_offset_no_lsb (struct pucr_node graph[], uint16_t in_len,
					       uint8_t no_esc,
					       uint8_t *no_lz_offset_lsb, uint8_t *no_lz2_offset_bits,
					       uint8_t *lz2_opt_huff_level,
					       uint8_t *lz_rank_prefix_len, uint8_t *lz_rank_prefix,
					       uint8_t *lz_rank_len, uint16_t lz_rank[]) {

	uint32_t lz_occ[16] = {0};
	uint16_t lz2_occ[16] = {0};

	//uint32_t lz_offset[65536]= {0};
	struct lz_offset_entry lz_offset[65536];
	for (uint i=0; i<65536; i++) {
		lz_offset[i].offset=i;
		lz_offset[i].count=0;
	}

	uint16_t lz2_full_occ[65536] = {0};

	// after all matches are finally determined, calculate the optimum of lsb-msb coding ratio
	// first step: count them, especially the bits their offsets need
	for (uint16_t i=0; i != in_len; i=graph[i].next_node) {
		if (graph[i].way_to_go == LZ) {
			// lz_length is not necessarily the real way to go,
			// but next_node is
			if ((graph[i].next_node - i) >= 3) {

				lz_offset[graph[i].lz_cur_offset].count++;

				// for the 3-byte and longer matches: calculate the case with "p" pure LSB and the rest E.Gamma-coded
			        for (int p=0; p <16; p++) // 0 = 1-bit offset, 1 = 2-bit offset, ...
					lz_occ[p] += (p+1) + int_log2( ((graph[i].lz_cur_offset-3)>>(p+1))+1 )*2+1;
//					lz_occ[p] += (p+1) + int_log2( ((graph[i].lz_cur_offset-(graph[i].next_node-i))>>(p+1))+1 )*2+1;
			} else {
				lz2_full_occ[graph[i].lz_cur_offset-2]++;
				int p = int_log2 (graph[i].lz_cur_offset - 2);
				p = p<=15?p:15;
				lz2_occ[p]++;
			}
		}
	}
	uint32_t min_sum = 0xFFFFFFFF;
	uint8_t min_no_lsb = 0;
	for (uint8_t m=15; m!=0xFF; m--) {
		if (lz_occ[m] < min_sum) {
			min_no_lsb = m;
			min_sum = lz_occ[m];
		}
	}
	*no_lz_offset_lsb = min_no_lsb + 1;

	uint8_t max_p = 16 - 2 - no_esc - 1;
	uint16_t min_lz2_bits = int_log2 (in_len);
	min_lz2_bits = (min_lz2_bits <= 8)?min_lz2_bits:8;
	int32_t max_saving = 0;
	for (uint8_t p=0; p<max_p; p++) {
		// count the bits saved for each LZ2 possible (i.e. if offset fits in p bits)
		int32_t gained = 0;
		for (uint8_t q=0; q<=p; q++) gained += (lz2_occ[q] * (16 - (p+1 + no_esc +2)));
		// lost bits due to unused LZ2 opportunities as offset length is longer than p bits
		int32_t lost = 0;
		for (uint8_t q=p+1; q<max_p; q++) lost += (lz2_occ[q] * (16 - (p+1 + no_esc +2)));
		// weigh against each other and store if minimum
		if ((gained-lost) > max_saving)  {
			max_saving = gained-lost;
			min_lz2_bits = p + 1;
		}
	}
	*no_lz2_offset_bits = min_lz2_bits;

// !!! should those be taken in account during the loop above??? !!!
uint8_t lit_cost[256]={8};
uint32_t lz2_huff_savings = huffiness (lz2_full_occ,0,(1<<(min_lz2_bits))-1, lz2_opt_huff_level,0, lit_cost);

	if (*lz_rank_prefix_len) {
		uint16_t min_prefix[min_no_lsb+1];
		uint16_t prefix_occ[min_no_lsb+1];

		if (min_no_lsb != -1) {
			 uint16_t prefix_sum;

			 for (uint8_t prefix_len=0; prefix_len <= min_no_lsb; prefix_len++) {
				prefix_occ[prefix_len] = 0xFFFF;
				min_prefix[prefix_len] = 0xFFFF;

				for (uint16_t prefix=0; prefix < (1<<prefix_len); prefix++) {
					uint16_t prefix_sum = 0;
					for (uint32_t i=3; i < in_len; i += (1<<(min_no_lsb+1)) )
						for (uint32_t j=prefix<<(min_no_lsb+1-prefix_len);(j<((prefix+1)<<(min_no_lsb+1-prefix_len))); j++)
							prefix_sum += lz_offset[i+j-3].count;

					if (prefix_sum < prefix_occ[prefix_len]) {
						min_prefix[prefix_len] = prefix;
						prefix_occ[prefix_len] = prefix_sum;
					}
				}
			}
		}

		// !!! couldn't this be done in parallel to the prefix_sum determination (requires min_no_lsb) ??? !!!
		for (uint16_t i=0;i < in_len;i++)
			// use the offset to be encoded, i.e. '-3'
			lz_offset[i].bits_used = (min_no_lsb+1) + int_log2( ((i-3)>>(min_no_lsb+1))+1 )*2+1;

		qsort( lz_offset, in_len, 3*sizeof(uint16_t), compare_lz_offset );

		// !!! couldn't this all be drawn into the lz_prefix loop above to be accounted for when chiosing min_no_lsb ??? !!!
		int32_t min_total = 0x7FFFFFFF, local_min, total;
		for (uint8_t prefix_len=0; prefix_len < min_no_lsb; prefix_len++) {
			local_min = 0x7FFFFFFF;
			for (uint8_t enum_bits=1; enum_bits < 8; enum_bits++) { // !!! upper bound 8 ??? !!!
			// costs vs. saved
				uint32_t cost = prefix_occ[prefix_len]*enum_bits  + ((1<<enum_bits)-1)*(min_no_lsb+5) + 6; // !!! +5 is some default for MSB, +6 (roughly) for header
				uint32_t save = 0;
				for (uint8_t i=0; i < (1<<enum_bits)-1; i++)
					// before the additional qsort and preceeding reduction of bits_used (see below)the following line was
					// save += lz_offset[i].count * (lz_offset[i].bits_used - prefix_len - enum_bits);
					save += lz_offset[i].count * (lz_offset[i].bits_used - enum_bits);
				total = cost - save;
				if (total < local_min)
					local_min = total;
				else
				    break; // if rising again, don't look further, it won't get better
				if ( (total < 0) && (total < min_total) ) {
					min_total = total;
					*lz_rank_prefix_len = prefix_len;
					*lz_rank_prefix = min_prefix[prefix_len];
					*lz_rank_len = (1<<enum_bits)-1;
					for (uint8_t i=0; i < (1<<enum_bits)-1; i++)
						lz_rank[i] = lz_offset[i].offset;
				}
			 }
			for (uint8_t i=0;i < 255; i++)
				lz_offset[i].bits_used = (lz_offset[i].bits_used<=1)?0:(lz_offset[i].bits_used-1);
			qsort( lz_offset, 255, 3*sizeof(uint16_t), compare_lz_offset );
		}
		if (min_total > 0)
			*lz_rank_len = 0;
	}
	return (0);
}


// ---=== ESC ESC ESC ESC ESC ESC ESC ESC ESC ESC ESC ESC ESC ESC ESC ===---


/*
    The algorithm in the OptimizeEscape() works as follows:
    1) Only unpacked bytes are processed, they are marked
       with MMARK. We proceed from the end to the beginning.
       Variable A (old/new length) is updated.
    2) At each unpacked byte, one and only one possible
       escape matches. A new escape code must be selected
       for this case. The optimal selection is the one which
       provides the shortest number of escapes to the end
       of the file,
        i.e. A[esc] = 1+min(A[0], A[1], .. A[states-1]).
       For other states A[esc] = A[esc];
       If we change escape in this byte, the new escape is
       the one with the smallest value in A.
    3) The starting escape is selected from the possibilities
       and mode 0 is restored to all mode 3 locations.
 */

static uint32_t optimize_esc (struct pucr_node graph[], const uint8_t *inbuf, uint16_t in_len,
                              uint8_t no_esc, uint8_t *startEscape,
			      uint8_t mtf, uint8_t mtf_second_line) {

    uint16_t i;
    int16_t j;

    uint8_t esc8 = 8-no_esc;

    uint16_t states = (1<<no_esc);
    int32_t minp = 0, minv = 0, other = 0;
    int32_t a[256]; /* needs int/long */
    int32_t b[256]; /* Remembers the # of escaped for each */


    for (i=0; i<256; i++)
        b[i] = a[i] = -1;

    uint8_t alphabet[256]; for (i=0; i<256;i++) alphabet[i]=i;

    for (i=0; i != in_len; i = graph[i].next_node)
        if (graph[i].way_to_go == LIT) {
            graph[i].way_to_go = MARKED;
            // store the (maybe to be mtf-encoded) literal in the graph to preserve inbuf
            graph[i].cur_lit = move_to_front_encode (inbuf[i], alphabet, mtf_second_line);
        }

    for (i=in_len-1; i != 0xFFFF; i--) {
        if (graph[i].way_to_go == MARKED) {
            graph[i].way_to_go = LIT;
            // use mtf-encoded literal for ESCape code evaluation
            int16_t k = graph[i].cur_lit >> esc8;

            /*
                k are the matching bytes,
                minv is the minimum value of neccessary following escapes
                minp is the minimum index, i.e. the esc code that currently delivers the least escapes
             */
           graph[i].new_esc = (minp << esc8);

            /* with k as esc code, there would be needed one more escape than the so-far minimum.
               not for the others, though. */
            a[k] = minv + 1;
            /* escapes happened so far. as long as k != minp, b[k] remains stuck at b[minp]+1.
               if we see a hit (k == minp), counting continues in the new b[k] track (with minp := k) */
            b[k] = b[minp] + 1;

            if (k==minp) {
                /* this is a hit, it requires an escape */
                minv++;
                /* Minimum changed -> need to find a new minimum */
                /* a[k] may still be the minimum */
                for (k=states-1; k>=0; k--) {
                    if (a[k] < minv) {
                        minv = a[k];
                        minp = k;
                        /*
                            There may be others, but the first one that
                            is smaller than the old minimum is equal to
                            any other new minimum (because minv was
                            incremented by 1 only).
                         */
                        break;
                    }
                }
            }
        }
    }
   /* Select the best value for the initial escape */
    /* find smallest a[j] with smallest "largest" escape code value */
    if (startEscape) {
        i = in_len;      /* make it big enough */
        for (j=0; j < states; j++) {
            if (a[j] <= i) {
                *startEscape = (j << esc8);
                i = a[j];
            }
        }
    }

    return b[startEscape ? (*startEscape>>esc8) : 0] + 1;
}


// ---=== PATH FINDER   PATH FINDER   PATH FINDER   PATH FINDER   PATH FINDER ===---


static uint32_t find_best_path (struct pucr_node graph[], const uint8_t *inbuf, uint16_t in_len,
				uint16_t rle_occ[], uint8_t rle_rank[], uint8_t rle_rank_len,
				uint8_t no_esc,
				uint8_t lz_lsb, uint8_t lz2_bits,
				uint8_t lit_cost_arr[],
				uint8_t fast_lane) {



// !!! not for empty input -- sure ???
	/* create rle cost (per character) table w/o cost for following length value */
	uint16_t rle_cost_table[256];
	for (uint16_t j=0; j < 256; j++) {
		/* is it a ranked character? */
		uint8_t k;
		for (k=0; k < rle_rank_len; k++) if ( j == rle_rank[k] ) break;
		if (k != rle_rank_len) {
			// ranked
			rle_cost_table[j] = 3 + int_log2(k+1)*2+1; // + (32 / (3 * rle_occ[j] + 1)); // this is a try to take table costs into account !!!
		} else {
			// not ranked
			rle_cost_table[j] = 3 + int_log2(rle_rank_len)+1 + 8;
			if (rle_rank_len == 0) rle_cost_table[j]--;
		}
	}
	/* however, it is some blur left here: depending on choices in the following optimizer loop, rle occurences and their
	   cost may change again, which could affect optimizer choices again, which might change costs again, which might ... */

	// graph[in_len] is the (virtual) EOF symbol
	graph[in_len].bits_to_end = 0;
	/* check for each node the best way to go */
	for (uint16_t i=in_len-1; i!=0xFFFF ; i--) {
		uint32_t rle_cost, lz_cost, dict_cost, lit_cost;
		uint32_t min = 0xFFFFFFFF;
		uint16_t rle_way, lz_way, dict_way;

		/* determine possible rle costs */
		/* this loops also checks for all shorter lengths */
		/* !!! ref org line 1044: only shorter ones that are a power of two? */
		for (uint16_t j=graph[i].rle_count; j > 1; j=fast_lane?1<<(int_log2(j)-1):j-1 ) {
			rle_cost = rle_cost_table[inbuf[i]] + int_log2(j-1)*2+1;
			rle_cost += graph[i+j].bits_to_end;
			if (rle_cost < min) {
				min = rle_cost;
				rle_way = i+j;
			}
		}
		// only if a new minimum was found, add no_esc
		// otherwise, keep up the high cost of 0xFFFFFFFF
		// !!! (min+1) == 0 is a bit hacky and relies on data type properties
		rle_cost = (min+1)?min+no_esc:min;

		/* determine possible LZ cost */
		min = 0xFFFFFFFF;
		/* !!! ref org line 1044: only shorters that are a power of two? */
		for (uint16_t j=graph[i].lz_count; j > 1; j=fast_lane?1<<(int_log2(j)-1):j-1) {
			// in case j == 2, the offset to be encoded may not be larger than 2^lz2_bits; otherwise lz is no option
			if ( (j == 2) && ((graph[i].lz_max_offset-2) >= (0x01 << lz2_bits)) ) break;
			lz_cost = (j==2)?2+lz2_bits:
//				 lz_lsb + int_log2(((graph[i].lz_max_offset-3)>>lz_lsb)+1)*2+1 + int_log2( graph[i-graph[i].lz_cur_offset].next_node - (i-graph[i].lz_cur_offset) +1 -1 +1)*2+1;
				  int_log2(j-1)*2+1 + lz_lsb + int_log2(((graph[i].lz_max_offset-3)>>lz_lsb)+1)*2+1;
////				  int_log2(j-1)*2+1 + lz_lsb + int_log2(((graph[i].lz_max_offset-j)>>lz_lsb)+1)*2+1;
			lz_cost += graph[i+j].bits_to_end;
			if (lz_cost < min) {
				min = lz_cost;
				lz_way = i+j;
			}
		}
		// only if a minimum was found in the loop, add the number of ESCape bits
		// otherwise, keep up the high cost of 0xFFFFFFFF
		// !!! (min+1) != 0 is a bit hacky and relies on data type properties
		if (min+1) {
			lz_cost = min + no_esc;
			graph[i].lz_cur_offset = graph[i].lz_max_offset;
			// if using a shorter-than-max match, we might find a closer one
			// this will keep the output of the encoded current offset smaller
			// and also might lead to a lower number of bits used for encoding LZ2/LZ offsets
			if (fast_lane == 0) { // only in fast_lane 0
				uint16_t used_match_length = lz_way-i;
				if (used_match_length<graph[i].lz_count) {
					find_closer_match (graph, inbuf, in_len, i, used_match_length,i-graph[i].lz_max_offset+1);
					// recalculate costs
					lz_cost = (used_match_length==2)?2+lz2_bits:
						  int_log2(used_match_length-1)*2+1 + lz_lsb + int_log2(((graph[i].lz_cur_offset-3)>>lz_lsb)+1)*2+1;
//						  int_log2(used_match_length-1)*2+1 + lz_lsb + int_log2(((graph[i].lz_cur_offset-used_match_length)>>lz_lsb)+1)*2+1;
					lz_cost += graph[i+used_match_length].bits_to_end+no_esc;
				}
			}
		} else {
			lz_cost = min;
		}

//		lit_cost = lit_cost_arr[graph[i].cur_lit] +graph[i+1].bits_to_end ; // + 1 ??? !!! ;
//		lit_cost = 7 +graph[i+1].bits_to_end ; // 7 to repect entropy encoding ??? !!! ;
		lit_cost = 8 +graph[i+1].bits_to_end ; // + 1 ??? !!! ;

		if ( (rle_cost < lz_cost) && (rle_cost < lit_cost) ) {
			graph[i].way_to_go = RLE;
			graph[i].next_node = rle_way;
			graph[i].bits_to_end = rle_cost;
		} else if ( (lz_cost <= rle_cost) && (lz_cost < lit_cost) ) {
			graph[i].way_to_go = LZ;
			graph[i].next_node = lz_way;
			graph[i].bits_to_end = lz_cost;
		} else {
			graph[i].way_to_go = LIT;
			graph[i].next_node = i+1;
			graph[i].bits_to_end = lit_cost;
		}
	}
fprintf (stderr, "--- PATH FINDER BITS TO END: %u \n",graph[0].bits_to_end);
}


// ---=== OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT ===---


static uint16_t outPointer = 0;
static uint8_t bitMask = 0x80;
static uint32_t bitCount = 0;

static int8_t out_buffer [OUT_SIZE];


static void flush_bits (void) {
    if (bitMask != 0x80) outPointer++;
}


static void put_bit (int bit) {
bitCount++;

    if (bit && outPointer < OUT_SIZE)
        out_buffer[outPointer] |= bitMask;
    bitMask >>= 1;
    if (!bitMask) {
        bitMask = 0x80;
        outPointer++;
    }
}


static void put_value (int value, int maxGamma) {
    int bits = 0, count = 0;

// !!! not working properly for maxGamma == 00 (output '0' or '1') !!! ???

    while (value>1) {
        bits = (bits<<1) | (value & 1); /* is reversed compared to value */
        value >>= 1;
        count++;
        put_bit (0);
    }
    if (count<maxGamma)
        put_bit (1);
    while (count--) {
        put_bit ((bits & 1));     /* output is reversed again -> same as value */
        bits >>= 1;
    }
}


static void put_n_bits (int byte, int bits) {
    while (bits--)
        put_bit ((byte & (1<<bits)));
}


static void put_huffed_value (int value, int bits, uint8_t rec) {

	if (rec == 0) {
		put_n_bits(value, bits);
		return;
	}

	uint16_t mask = 1<<(bits-1);
	mask |= mask >> 1;

	if ((value & mask) == 0) {
		put_bit(0);
		put_huffed_value (value, bits-2, rec-1);
		return;
	}
	else if ((value & mask) == (01<<bits-2)) {
		put_bit (1);
		put_bit (0);
		put_n_bits (value, bits-2);

	} else {
		put_bit (1);
		put_bit (1);
		put_n_bits (value, bits-1);
	}
}

static void put_huffed_value_org (int value, int bits, uint8_t rec) {

	if (rec == 0) {
		put_n_bits(value, bits);
		return;
	}

	uint16_t mask = 1<<(bits-1);
	mask |= mask >> 1;

	if ((value & mask) == 0)
		if (rec == 1) {

		} else {
			put_bit(0);
			put_huffed_value (value, bits-2, rec-1);
			return;
		}
	else if ((value & mask) == (1<<bits-2)) {
		put_bit (1);
	} else {
		put_bit (1);
		put_bit (0);
	}
	put_n_bits (value, bits-1);
}

static uint32_t output_path (struct pucr_node graph[], const uint8_t *inbuf, uint16_t in_len,
			uint16_t *out_len,
			uint16_t rle_occ[], uint8_t rle_rank[], uint8_t rle_rank_len,
			uint8_t start_esc, uint8_t no_esc, uint8_t mtf, uint8_t mtf_second_line, uint8_t lit_opt_huff_level,
			uint8_t lz_lsb, uint8_t lz2_bits, uint8_t lz_length_max_gamma, uint8_t lz2_opt_huff_level,
			uint8_t lz_rank_prefix_len, uint8_t lz_rank_prefix, uint8_t lz_rank_len, uint16_t lz_rank[]) {


	uint8_t graph_current_esc = start_esc;
	uint8_t current_esc = start_esc;
	uint8_t current_esc_mask = 0xFF00 >> no_esc;
	uint8_t current_esc_8 = current_esc >> (8 - no_esc);

	// output header
	// decompressed length
	put_value ((in_len>>8)+1, NO_MAX_GAMMA);
	put_n_bits (in_len,8);
	// lz length max gamma
	put_n_bits (lz_length_max_gamma,4);
	// number of ESCape bits and first ESCape sequence
	put_n_bits (no_esc,4);
	put_n_bits (current_esc_8,no_esc);
	// RLE table count and table w/ up to 15 entries
	put_value (rle_rank_len+1, 6);
	for (uint8_t i=0; i < rle_rank_len; i++)
		put_n_bits(rle_rank[i],8);
	// number of bits for LZ2 offset
	put_n_bits (lz2_bits,4);
	// number of LSBs for LZ offset
	put_n_bits (lz_lsb,4);
	// LZ offset ranking
	uint8_t lz_rank_len_bits = int_log2 (lz_rank_len+1);
	put_value (lz_rank_len_bits+1, NO_MAX_GAMMA);
	if (lz_rank_len > 0) {
		put_n_bits (lz_rank_prefix_len,4);
		put_n_bits (lz_rank_prefix, lz_rank_prefix_len);
		for (uint8_t i=0; i< lz_rank_len; i++) {
			put_n_bits (lz_rank[i]-3, lz_lsb);
			put_value (((lz_rank[i]-3) >> lz_lsb)+1, NO_MAX_GAMMA);
		}
	}
	// move to front
	put_value (mtf_second_line+1,NO_MAX_GAMMA);
	// re-think sharply: the number of bits to be encoded is limited by the number of esc-bits (and thus remaining lit-bits) !!!
	put_n_bits (lit_opt_huff_level, 3); // could reach up to 4 on case of 0 esc bits
	// re-think sharply: the number of bits to be encoded is limited by the number of lz2 offset-bits !!!
	put_n_bits (lz2_opt_huff_level,4);
	// derive further parameters
	uint8_t no_rle_char_esc_bits = (rle_rank_len==0)?0 : int_log2(rle_rank_len)+1;

	// output data
	for (uint16_t i=0; i!=in_len; i=graph[i].next_node) {
		switch (graph[i].way_to_go) {
			case RLE:
				// ESC + 100
				put_n_bits (current_esc_8, no_esc);
				put_bit (1); put_bit (0); put_bit (0);
                                // output RLE count (real count minus 1, i.e. 1 instead of 2 etc.) E. Gamma coded
                                put_value (graph[i].next_node - i - 1, NO_MAX_GAMMA);
				// check if it is a ranked character
				uint8_t k;
				for (k=0; k < rle_rank_len; k++) if ( inbuf[i] == rle_rank[k] ) break;
				if (k != rle_rank_len) {
					// is ranked: E. Gamma coded pointer to rank list
					put_value (k+1, NO_MAX_GAMMA);
				} else {
					// not ranked: determine '00..00'-code length (already determined outside loop) followed by character
					put_n_bits (0x00, no_rle_char_esc_bits);
					put_n_bits (inbuf[i], 8);
				}
				break;
			case LZ:
				put_n_bits (current_esc_8, no_esc);

uint16_t lz= graph[i].lz_cur_offset-(graph[i].next_node -i);

				if ( (graph[i].next_node-i) == 2 ) {
					// LZ2: ESC + 11 + LZ2 OFFSET (lz2_bits)
	                                put_bit (1); put_bit (1); // put_bit (1);
					if (lz2_bits <= 2)
						put_n_bits (graph[i].lz_cur_offset-2, lz2_bits);
					else
						put_huffed_value (graph[i].lz_cur_offset-2, lz2_bits,lz2_opt_huff_level);
				} else {
					// LZ: ESC + E. Gamma coded actually used lz length (next_node - i - 1), not neccessariliy lz_max_len ...
					put_value (graph[i].next_node - i - 1 , lz_length_max_gamma);
					// ... and also E. Gamma coded LZ offset minus the actually used length of lz match(=next_node-i)
					uint16_t off = (graph[i].lz_cur_offset);
					// do we need to take care of ranked offsets?
					uint8_t l;
					for (l=0; l < lz_rank_len; l++)
						if (off == lz_rank[l]) break;
					if (l != lz_rank_len) {
						// break, i.e. ranked offset
						put_n_bits (lz_rank_prefix, lz_rank_prefix_len);
						put_n_bits (l, lz_rank_len_bits);
					} else {
						// unranked offset
						off -= 3;
						// LSBs
						if ( (lz_rank_len) && (off & (0xFFFF >> (16-lz_lsb))) >> (lz_lsb - lz_rank_prefix_len) == lz_rank_prefix ) {
							put_n_bits (lz_rank_prefix, lz_rank_prefix_len);
							put_n_bits (0xFFFF, lz_rank_len_bits);
							put_n_bits (off, lz_lsb - lz_rank_prefix_len);
						} else {
							put_n_bits (off, lz_lsb);
						}
						// MSBs (+1 for E. Gamma code)
						put_value ((off >> lz_lsb)+1, NO_MAX_GAMMA);
					}
				}
				break;
			case LIT:
				if (graph[i].new_esc != graph_current_esc) {
//					fprintf (stderr, "  NEW possible ESC: %u", graph[i].new_esc);
					graph_current_esc = graph[i].new_esc;
				}
				if ( (graph[i].cur_lit & current_esc_mask) != current_esc) {
					// output the huffed mtf-encoded literal
					put_n_bits (graph[i].cur_lit >> (8-no_esc), no_esc);
					put_huffed_value(graph[i].cur_lit, (8-no_esc),lit_opt_huff_level);
				} else {
					// ESC = first bits of LITeral
					put_n_bits (current_esc_8, no_esc);
					// 101
	                                put_bit (1); put_bit (0); put_bit (1);
					// new ESC code
					current_esc = graph[i].new_esc;
					current_esc_8 = current_esc >> (8 - no_esc);
					put_n_bits (current_esc_8, no_esc);
					// rest of mtf-encoded LITeral:
					put_huffed_value(graph[i].cur_lit, (8-no_esc),lit_opt_huff_level);
				}
				break;
			default:
				// !!! shall not happen
				break;
		}
	}
	flush_bits();
	*out_len = (bitCount+7) >> 3;
}


void print_graph_stats (struct pucr_node graph[], uint8_t inbuf[], uint16_t in_len, uint8_t no_esc, uint8_t start_esc) {

	uint16_t i = 0;

	uint16_t rle_cnt = 0;
	uint16_t rle_bytes = 0;

	uint16_t lit_cnt = 0;
	uint16_t lit_esc_cnt = 0;
	uint16_t lit_unesc_cnt = 0;
	uint16_t cur_lit_occ[256] = {0};

	uint16_t unesc_seq_len = 0;
	uint16_t unesc_seq_cnt = 0;
	uint16_t unesc_seq_lit = 0;

	uint16_t esc_seq_len = 0;
	uint16_t esc_seq_cnt = 0;
	uint16_t esc_seq_x = 0;

	uint16_t lz_cnt = 0;
	uint16_t lz2_cnt = 0;
	uint16_t lz3_cnt = 0;
	uint16_t lz4_cnt = 0;
	uint16_t lz5_cnt = 0;
	uint16_t lz6_cnt = 0;
	uint16_t lz7_cnt = 0;
	uint16_t lz8_cnt = 0;
	uint16_t lz_bytes = 0;

	uint16_t twogra_cnt;

	uint16_t lz_first_match = 0;

	uint16_t non_ip_cnt = 0;
	uint16_t ip_cnt = 0;
	uint32_t sum_ip = 0;

	uint16_t ip_corr_cnt = 0;

	uint16_t lz_off_cnt = 0;
	uint32_t sum_lz_off = 0;

	uint8_t current_esc = start_esc;
	uint8_t graph_current_esc = start_esc;

	uint8_t current_esc_mask = 0xFF00 >> no_esc;

	uint8_t min_esc_seq_len = 16 >> no_esc;

	uint16_t n = 0;
	uint16_t ip[65536] = {0xFFFF};
	for (i=0; i != in_len; i=graph[i].next_node) {
		uint8_t esc = 1;
		uint8_t lit = 0;

// if (graph[i].way_to_go != RLE)
ip[n++] = i;


		if ( graph[i].way_to_go == RLE) {
			rle_cnt++;
			rle_bytes += graph[i].next_node - i;
		}
		if ( graph[i].way_to_go == LZ) {
			lz_cnt++;
			uint16_t dist = graph[i].next_node - i;
			if (dist == 2) lz2_cnt++;
			if (dist == 3) lz3_cnt++;
			if (dist == 4) lz4_cnt++;
			if (dist == 5) lz5_cnt++;
			if (dist == 6) lz6_cnt++;
			if (dist == 7) lz7_cnt++;
			if (dist == 8) lz8_cnt++;
			lz_bytes += dist;
			uint16_t off = i-graph[i].lz_cur_offset;
			if (dist > 2 && n > 2) {
				// is offset in array of interesting points?
				uint16_t j;
				for (j=0; j<(n-1); j++) if (off == ip[j]) break;
				if (j == (n-1) ) { // not in array
					uint16_t jj;
					for (jj = (n-2); jj > 0; jj--) {
if (ip[jj]+dist>=i) continue; // no overlap
						uint16_t k;
						for (k = 1; k < dist; k++) if (inbuf[i+k] != inbuf[ip[jj]+k]) break;
						if (k == dist) {
							ip_corr_cnt++;
							ip_cnt++;
							sum_ip += (n-jj-1);
// graph[i].way_to_go = LZ_IP;
// graph[i].lz_cur_offset=n+3-jj; // i-ip[jj];
							// sum_lz_off += graph[i].lz_cur_offset-3;
							// lz_off_cnt++;
							break;
						}
					}
					if (jj == 0) {
						non_ip_cnt++;
						sum_lz_off += graph[i].lz_cur_offset-3;
						lz_off_cnt++;
					}
				} else {
					ip_cnt++;
					sum_ip += (n-j-1);
// graph[i].way_to_go = LZ_IP;
// graph[i].lz_cur_offset=n+3-j; // i-ip[j]; // unneccessary, remains same
				}
				if ( graph[i].seen_before == off ) lz_first_match++;

			}
		}
		if ( graph[i].way_to_go == LIT ) {
			lit = 1;
			lit_cnt++;
			cur_lit_occ[graph[i].cur_lit]++;
			if (graph[i].new_esc != graph_current_esc) {
				// new possible (!) ESC code if need to switch
				graph_current_esc = graph[i].new_esc;
			}
			if ( (graph[i].cur_lit & current_esc_mask) != current_esc) {
				esc = 0;
			} else {
				lit_esc_cnt++;
				current_esc = graph[i].new_esc;
			}
		} else {
			// any non-LIT
		}

		if (!esc) {
			unesc_seq_len++;
			if (esc_seq_len > min_esc_seq_len) {
				esc_seq_cnt++;
				esc_seq_x += esc_seq_len;
			}
			esc_seq_len = 0;
		} else {
			esc_seq_len++;
			if (unesc_seq_len > 8) {
				unesc_seq_cnt++;
				unesc_seq_lit += unesc_seq_len;
			}
			unesc_seq_len = 0;
		}
	}
	lit_unesc_cnt = lit_cnt - lit_esc_cnt;

fprintf(stderr, "\
------== GRAPH  STATS ==------\n\
#RLE CNT: 		%u\n\
#RLEed BYTES:		%u\n\
#LIT:			%u\n\
 #ESC:		%u\n\
 #UNESC:	%u\n\
#UNESC LIT SEQs >8:	%u\n\
#UNESC LITs IN SEQs	%u\n\
#ESC SEQs >%u:		%u\n\
#ESCs IN SEQs		%u\n\
#LZ  CNT:		%u\n\
#LZ2 CNT:		%u\n\
#LZ3 CNT:		%u\n\
#LZ4 CNT:		%u\n\
#LZ5 CNT:		%u\n\
#LZ6 CNT:		%u\n\
#LZ7 CNT:		%u\n\
#LZ8 CNT:		%u\n\
LZed BYTES:             %u\n\
LZ FIRST MATCH		%u\n\
#IP ENTRIES		%u\n\
NON IP REFS		%u\n\
IP REF HITS		%u\n\
IP CORR CNT		%u\n\
IP AVG INDX		%u\n\
LZ OFF AVG		%u\n\
TWOGRAPHs:		%u\n\
------------------------------\n",
rle_cnt, rle_bytes, lit_cnt,lit_esc_cnt,lit_unesc_cnt, unesc_seq_cnt,unesc_seq_lit, min_esc_seq_len,esc_seq_cnt,esc_seq_x, lz_cnt,lz2_cnt,lz3_cnt,lz4_cnt,lz5_cnt,lz6_cnt,lz7_cnt,lz8_cnt,lz_bytes,lz_first_match,n,non_ip_cnt,ip_cnt,ip_corr_cnt,sum_ip/ip_cnt,sum_lz_off/lz_off_cnt,twogra_cnt);

//for (i=0; i< 256; i++) {
// fprintf (stderr, "LIT %3u: %5u x\n",i, cur_lit_occ[i]);
//}

// couldn't this be returned by update_lit_occ_from_graph?!
uint16_t sum = 0;
for (int16_t i = 0; i< 256; i++) sum += cur_lit_occ[i];


}


uint32_t pucrunch_256_encode (uint8_t *outbuf,uint16_t *out_len,
                             uint8_t *inbuf, uint16_t in_len, // const !!! inbuf
			     uint8_t fast_lane) {

start_clock();

	struct pucr_node graph[65536];
	uint16_t rle_occ[256];
	uint8_t rle_rank[MAX_RLE_RANKS];
	uint8_t rle_rank_len;

	uint16_t lit_occ[256];
	uint8_t lit_cost[256] = {8};
	uint8_t opt_level;

	uint8_t lit_rank[256];
	uint8_t lit_rank_len;

	count_rle_and_generate_occ_from_inbuf (graph, inbuf, in_len, rle_occ);
	generate_rle_ranks (rle_occ, rle_rank, &rle_rank_len);
fprintf (stderr, "--- COUNTING RLE: ");
print_clock(); start_clock();

uint8_t lz_length_max_gamma = LZ_LEN_MAX_GAMMA;

	uint8_t lz_lsb = int_log2 (in_len)*2/3;
	uint8_t lz2_bits = int_log2 (in_len)*2/3;
		lz2_bits = lz2_bits>8?8:lz2_bits;
	uint8_t lz2_opt_huff_level = 0;

	uint32_t min_no_bits = 0xFFFFFFFF;
	uint8_t min_no_esc;
	uint8_t min_lz_lsb;
	uint8_t min_lz2_bits;
	uint8_t min_lz2_opt_huff_level = 0;

	uint8_t min_start_esc;

        prepare_last_seen_from_graph (graph, inbuf, in_len);
	find_matches (graph, inbuf, in_len, 0, in_len-1, 2, fast_lane, lz_length_max_gamma, lz2_bits, lz_lsb);
//	find_dict (graph, inbuf, in_len);

uint8_t lz_rank_prefix_len, lz_rank_prefix, lz_rank_len;
uint8_t min_lz_rank_prefix_len, min_lz_rank_prefix, min_lz_rank_len;
uint16_t lz_rank[127];
uint16_t min_lz_rank[127];

	// do the graph building with (standard) mtf turned on
	uint8_t mtf = 1;
if (fast_lane > 1) mtf = 0;


	int8_t no_esc;
	uint8_t start_esc = 0;

	if (fast_lane < 3) {
	for (no_esc = (fast_lane<2)?0:1; no_esc < ((fast_lane<2)?MAX_ESC:4); no_esc++) {
			// the first path generated using sane default

lz_lsb = int_log2 (in_len)*2/3;
lz2_bits = int_log2 (in_len)*2/3;
lz2_bits = lz2_bits>8?8:lz2_bits;
for (int i=0;i<256;i++) lit_cost[i]=8;
			find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rle_rank_len, no_esc, lz_lsb, lz2_bits, lit_cost, fast_lane);
// fprintf (stderr, "AFTER FIRST OPTIMIZE: %u bits\n",graph[0].bits_to_end);
			lz_rank_prefix_len = 0;
			find_optimal_lz_offset_no_lsb (graph, in_len, no_esc, &lz_lsb, &lz2_bits, &lz2_opt_huff_level, &lz_rank_prefix_len, &lz_rank_prefix, &lz_rank_len,lz_rank);

// fprintf (stderr, "NO ESC: %u --- LZ2: %u --- LZLSB: %u \n",no_esc,lz2_bits,lz_lsb);
			int escaped = optimize_esc (graph, inbuf, in_len, no_esc, &start_esc, mtf, 7);

if (fast_lane < 2)	update_rle_occ_from_graph (graph, inbuf, in_len, rle_occ);
if (fast_lane < 2)	generate_rle_ranks (rle_occ, rle_rank, &rle_rank_len);

if (!fast_lane)		update_lit_occ_from_graph (graph, in_len, lit_occ);
if (!fast_lane)		huffiness (lit_occ, 0, 0xFF>>no_esc, &opt_level,1,lit_cost);

if (fast_lane < 2)	find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rle_rank_len, no_esc, lz_lsb, lz2_bits, lit_cost, fast_lane);

// fprintf (stderr, "AFTER SECOND OPTIMIZE: %u bits\n",graph[0].bits_to_end);
if (fast_lane < 2) {	lz_rank_prefix_len = 1;
			find_optimal_lz_offset_no_lsb (graph, in_len, no_esc, &lz_lsb, &lz2_bits, &lz2_opt_huff_level, &lz_rank_prefix_len, &lz_rank_prefix, &lz_rank_len,lz_rank);
		   }
// find_optimal_lz_offset_no_lsb (graph, in_len, no_esc, &lz_lsb, &lz2_bits, &lz2_opt_huff_level,0,0,0,0);

	            	/* Compare value: bits lost for escaping -- bits lost for prefix */
			if ((graph[0].bits_to_end + (no_esc+3) * escaped) < min_no_bits) {
				min_no_esc = no_esc;
				min_no_bits = graph[0].bits_to_end + (no_esc+3)*escaped; // "- escaped" is for those counted 9-bit LIT in find_best_path
				// also store the coresponding lz_lsb value
				min_lz_lsb = lz_lsb;
				min_lz2_bits = lz2_bits;
				min_lz2_opt_huff_level = lz2_opt_huff_level;
				min_lz_rank_len = lz_rank_len;
				min_lz_rank_prefix_len = lz_rank_prefix_len;
				min_lz_rank_prefix = lz_rank_prefix;
				for (uint8_t i=0; i < lz_rank_len; i++)
					min_lz_rank[i] = lz_rank[i];
				min_start_esc = start_esc & ((0xff00>>no_esc) & 0xff);
			}
fprintf (stderr, "--- OPTIMIZING ESC %u: ", no_esc);
print_clock(); start_clock();
		}
//	    }
	} else {
		// fast_lane >= 3:
		// just assume some sane defaults from above for further use
		min_lz_lsb = lz_lsb;
		min_lz2_bits = lz2_bits;
		min_no_esc = 2;
		min_start_esc = 0;
		min_lz_rank_len = 0;
	}
for (int i=0;i< 256;i++) lit_cost[i]=8;
	// (re-)calculate best path with the values found or set above
	// alternatively, store the best graph found in the loop so far... (memory!!!)
	find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rle_rank_len, min_no_esc, min_lz_lsb, min_lz2_bits, lit_cost, fast_lane);
	int escaped = optimize_esc (graph, inbuf, in_len, min_no_esc, &start_esc, mtf, 7);

if (fast_lane < 2)		update_rle_occ_from_graph (graph, inbuf, in_len, rle_occ);
if (fast_lane < 2)		generate_rle_ranks (rle_occ, rle_rank, &rle_rank_len);

 if (!fast_lane)		update_lit_occ_from_graph (graph, in_len, lit_occ);
 if (!fast_lane)		huffiness (lit_occ, 0, 0xFF>>no_esc, &opt_level,1,lit_cost);




if (fast_lane < 2)		find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rle_rank_len, min_no_esc, min_lz_lsb, min_lz2_bits, lit_cost, fast_lane);
// if (fast_lane < 2)		update_rle_occ_from_graph (graph, inbuf, in_len, rle_occ);
// if (fast_lane < 2)		generate_rle_ranks (rle_occ, rle_rank, &rle_rank_len);

// fprintf (stderr, " LZ PREFIX LEN = %u  LZ PREFIX = %u  LZ RANKED = %u\n",min_lz_rank_prefix_len, min_lz_rank_prefix, min_lz_rank_len);

	// handling of remaining LITerals:
	int32_t huffed = 0;

fprintf (stderr, "GENERATING FINAL GRAPH: ");
print_clock(); start_clock();

	// if mtf, will a higher mtf_second_line be beneficial?

// min_lit_opt_huff_level = (8 - min_no_esc) >> 1;
// min_mtf_second_line = (8 - min_no_esc) >> 1;

uint8_t min_mtf = 1;
	uint8_t min_mtf_second_line = 0;
	uint8_t  min_lit_opt_huff_level = 0;
	uint32_t min_mtf_second_line_bits = 0xFFFFFFFF;
		for (uint8_t mtf_second_line=0; mtf_second_line < (fast_lane?2:MAX_MTF); mtf_second_line++) {
			escaped = optimize_esc (graph, inbuf, in_len, min_no_esc, &min_start_esc, min_mtf, mtf_second_line);
			min_no_bits = graph[0].bits_to_end + (min_no_esc+3)*escaped;

			huffed=0;
			update_lit_occ_from_graph (graph, in_len, lit_occ);
			huffed = huffiness (lit_occ, 0, 0xFF >> min_no_esc , &opt_level,1,lit_cost);

  fprintf (stderr, "MTF CHECK %u BITS: %u\n",mtf_second_line, min_no_bits-huffed);
  fprintf (stderr,"huffed = %u @ level %u\n",huffed,opt_level);

			if ( (min_no_bits-huffed) < min_mtf_second_line_bits ) {
				min_mtf_second_line = mtf_second_line;
				min_mtf_second_line_bits = (min_no_bits - huffed);
				min_lit_opt_huff_level = opt_level;
			}
		}

fprintf (stderr, "----- MTF: %u    MTF SECOND LINE: %u\n",min_mtf, min_mtf_second_line);
fprintf (stderr,"OPT LIT HUFF LEVEL: %u\n",min_lit_opt_huff_level);

fprintf (stderr, "HUFFMAN ON LITERALS: ");
print_clock(); start_clock();
	escaped = optimize_esc (graph, inbuf, in_len, min_no_esc, &min_start_esc, min_mtf, min_mtf_second_line);
huffed=0;

fprintf (stderr, "CALCULATE (AGAIN) USING FINALLY DETERMINED PARAMETERS: ");
print_clock(); start_clock();


fprintf (stderr, "START ESC: %u\n",min_start_esc);
fprintf(stderr, 
"-----== ENCODING STATS ==-----\n#ESC: 			%u\n#ESCAPED LIT:		%u\n#LZ OFFSET LSB:		%u\n\
#LZ2 OFFSET BITS	%u\n#RANKED RLEs		%u\nIN:			%u\nOUT: without header	%u\n     \
before huffman\n------------------------------\n",
min_no_esc, escaped, min_lz_lsb, min_lz2_bits, rle_rank_len, in_len, (min_no_bits + 7)>>3 );

// print_graph_stats (graph, inbuf, in_len, min_no_esc, min_start_esc);
// find_optimal_lz_offset_no_lsb (graph, in_len, no_esc, &lz_lsb, &lz2_bits, &lz2_opt_huff_level, &lz_rank_prefix_len, &lz_rank_prefix, &lz_rank_len,lz_rank);
// find_optimal_lz_offset_no_lsb (graph, in_len, min_no_esc, &min_lz_lsb, &min_lz2_bits, &min_lz2_opt_huff_level, &lz_rank_prefix_len, &lz_rank_prefix, &lz_rank_len,lz_rank);
// fprintf(stderr, "-----== AFTER ENCODING STATS ==-----\n#LZ OFFSET LSB:		%u\n#LZ2 OFFSET BITS	%u\n------------------------------\n",lz_lsb, lz2_bits);
// uint8_t lz_ai_offset = 0;
// find_optimal_lz_ip_offset_no_lsb (graph, in_len, &lz_ai_offset);
// min_lz_lsb--;

	output_path (graph, inbuf, in_len, out_len, rle_occ, rle_rank, rle_rank_len, 
			min_start_esc, min_no_esc, 
			min_mtf, min_mtf_second_line, min_lit_opt_huff_level, 
			min_lz_lsb, min_lz2_bits, lz_length_max_gamma, min_lz2_opt_huff_level, 
			min_lz_rank_prefix_len, min_lz_rank_prefix, min_lz_rank_len,min_lz_rank); 

fprintf (stderr, "--- OUTPUT FINAL PATH: ");
print_clock(); start_clock();
}

// ---=== DECODE DECODE DECODE DECODE DECODE ===---

// static const uint64_t *up_Data;
static const uint8_t *up_Data;
static uint64_t upd;
static uint64_t up_Mask;
// static const unsigned char *up_Data;
// static uint32_t up_Mask;
// static uint16_t up_Byte;
 static uint32_t up_bitCount;

static void up_SetInput(const void *data) {
// static void up_SetInput(const unsigned char *data) {
    up_Data = data;
upd = be64toh(*up_Data);
    up_Mask = 0x8000000000000000;
//    up_Mask = 0x80;
//    up_Byte = 0;
    up_bitCount = 0;
}


static inline uint32_t rotl32 (uint32_t n, uint32_t c)
{
  const unsigned int mask = (CHAR_BIT*sizeof(n) - 1);

  c &= mask;
  return (n<<c) | (n>>( (-c)&mask ));
}

static inline uint64_t rotl64 (uint64_t n, uint32_t c)
{
  const unsigned int mask = (CHAR_BIT*sizeof(n) - 1);

  c &= mask;
  return (n<<c) | (n>>( (-c)&mask ));
}

static inline uint64_t rotl64_alt ( uint64_t x, uint8_t r )
{
  return (x << r) | (x >> (64 - r));
}

static inline uint32_t rotl32_asm (uint32_t u, size_t r)
{
    __asm__ ("roll %%cl, %0" : "+r" (u) : "c" (r));
    return u;
}

const uint32_t mask[32] = {0,1,3,7,15,31,63,127,255,511,1023,2047,4095,8191,16383,32767,65535,131071,262143,524287,1048575,2097151,4194303,8388607};

uint32_t inp;

static uint32_t inline up_GetBits(uint32_t bits) {

// uint64_t inp = *(up_Data + (up_bitCount >> 3));
// uint64_t inp = up_Data[up_bitCount >> 3];

	uint32_t inp;

	//	memcpy (&inp, up_Data+(up_bitCount >> 3),4);
	uint32_t *buf;
	buf = &up_Data[up_bitCount >>3];
	inp = *buf;

	inp = be32toh(inp);
	inp = rotl32(inp, (up_bitCount & 7) + bits);
//	inp >>= 32- ((up_bitCount & 7) + bits);
	up_bitCount+=bits;
	return inp & ((1 << bits)-1);
//  return inp & mask[bits];
}

static uint32_t inline up_GetMoreBits(uint32_t bits) {

// uint64_t inp = *(up_Data + (up_bitCount >> 3));
// uint64_t inp = up_Data[up_bitCount >> 3];

/*	uint32_t inp;

	//	memcpy (&inp, up_Data+(up_bitCount >> 3),4);
	uint32_t *buf;
	buf = &up_Data[up_bitCount >>3];
	inp = *buf;

	inp = be32toh(inp);
	inp = rotl32(inp, (up_bitCount & 7) + bits);
//	inp >>= 32- ((up_bitCount & 7) + bits);
*/
	inp = rotl32(inp, bits);
	up_bitCount+=bits;
	return inp & ((1 << bits)-1);
//  return inp & mask[bits];
}

static void inline up_ConsumeBits (uint32_t bits) {
	up_bitCount += bits;
}


static uint32_t inline up_PeekBits(uint32_t bits) {

// uint64_t inp = *(up_Data + (up_bitCount >> 3));
// uint64_t inp = up_Data[up_bitCount >> 3];

	uint32_t inp;

	//	memcpy (&inp, up_Data+(up_bitCount >> 3),4);
	uint32_t *buf;
	buf = &up_Data[up_bitCount >>3];
	inp = *buf;

	inp = be32toh(inp);
	inp = rotl32(inp, (up_bitCount & 7) + bits);
	return inp & ((1 << bits)-1);
}


static uint32_t inline up_GetOneBit(void) {

	uint8_t inp;

	//	memcpy (&inp, up_Data+(up_bitCount >> 3),1);
	uint8_t *buf;
	buf = &up_Data[up_bitCount >>3];
	inp = *buf;

	inp &= (0x80 >> (up_bitCount & 7));
	up_bitCount++;
	return inp;
}


static uint32_t inline up_GetValue(void) {

// this could be faster in 32-bit version
// 64 bit is needed only for long files with
// long RLE (delenn.bin) or maybe LZ.
// not so much fpr shorter ethernet packets.
// consider 32-bit then or restrict value size

	uint64_t inp;

 	// memcpy (&inp, up_Data+(up_bitCount >> 3),4);
	uint64_t *buf;
	buf = &up_Data[up_bitCount >>3];
	inp = *buf;

	inp = be64toh(inp);
	inp <<= up_bitCount & 7;
	uint32_t n = __builtin_clzll(inp);

	uint32_t m = (n<<1) +1;
	inp = rotl64(inp, m);
	up_bitCount += m;

	// shifting by (n+1) would be sufficent; instead,
	// m is used as it already is available
	return (inp & ((1 << m)-1));
}


static uint32_t inline up_GetValueOrg(void) {
    uint8_t i = 0;

    while (1) {
//        if (unlikely(!up_GetBits(1)))
        if (unlikely(up_GetOneBit()))
            break;
        i++;
    }

    return (1<<i) | up_GetBits(i);
}


inline unsigned int gt(uint32_t x, uint32_t y) {
	return (y-x) >> (CHAR_BIT*sizeof(int)-1);
}


static uint32_t inline up_GetValueMaxGamma (uint32_t maxGamma) {

	uint32_t inp;

	//	memcpy (&inp, up_Data+(up_bitCount >> 3),4);
	uint32_t *buf;
	buf = &up_Data[up_bitCount >>3];
	inp = *buf;

	inp = be32toh(inp);
	inp <<= up_bitCount & 7;

	uint32_t n = __builtin_clz(inp);
	n = min (n, maxGamma);
	uint32_t m = gt (maxGamma,n);
	inp = rotl32(inp, 2*n+1);

	inp >>=1-m;
	up_bitCount += (n<<1)+m;

	return (1<<n) | (inp & ((1 << n) -1));
}


static uint32_t inline up_GetValueMaxGammaOrg(uint32_t maxGamma) {
    uint8_t i = 0;

    while (i<maxGamma) {
        if (unlikely(up_GetBits(1)))
//        if (unlikely(!up_GetOneBit()))
            break;
        i++;
    }
    return (1<<i) | up_GetBits(i);
}


static uint32_t inline up_get_huffed_value (uint32_t bits, int32_t lvl) {

	uint32_t inp;

	//	memcpy (&inp, up_Data+(up_bitCount >> 3),4);
	uint32_t *buf;
	buf = &up_Data[up_bitCount >>3];
	inp = *buf;

	inp = be32toh(inp);
	inp <<= up_bitCount & 7;

	int32_t n = __builtin_clz(inp);
	n = min (n, lvl);

	uint32_t m = gt (lvl,n); // '1' if lvl != n, '0' if lvl == n

	uint32_t r = (bits - (n<<1));
	r -= (m << 1); // -0 if lvl or -2 otherwise

	inp = rotl32(inp, n+(m<<1)); // n plus 0 or 2 bits more

	r += (m & inp); //  '0' for left sub-tree (10) or '1' for right sub-tree (11)

	up_bitCount += n+(m<<1)+r;

	inp = rotl32(inp,r);
	inp &= (1 << r) - 1;
	inp |= m << r;
	return (inp);
}


static uint32_t inline up_get_lz_offset (uint32_t prefix_len, uint32_t prefix, uint32_t index_len, uint32_t offset_len, uint16_t ranked_lz_offset[]) {

	uint32_t inp;

//	memcpy (&inp, up_Data+(up_bitCount >> 3),4);
	uint32_t *buf;
	buf = &up_Data[up_bitCount >>3];
	inp = *buf;

	inp = be32toh(inp);
	inp <<= up_bitCount & 7;

	inp = rotl32(inp,prefix_len);
	uint32_t off = inp & ((1<<prefix_len)-1);

	uint32_t d = off ^ prefix;
		 d = sign_extension (d-1); // prefix ? --> d = '1...1'

	inp = rotl32(inp, index_len & d);
	uint32_t i = inp & ((1 << index_len)-1);
	up_bitCount += index_len & d;

        d &= sign_extension ( i- ((1 << index_len)-1) ); // real index? --> d = '1..1'

	if (likely(d==0)) { // okay, this is not branchless anymore...
			    // but completely branchless version was slower
		inp = rotl32(inp, offset_len - prefix_len);
		off <<= (offset_len - prefix_len);
		off += inp & ((1 << (offset_len - prefix_len))-1);

		uint32_t n = __builtin_clz(inp);
		uint32_t m = (n<<1) +1;
		inp = rotl32(inp, m);

		up_bitCount += ((offset_len + m) );

		off |= ( (inp & ((1 << (n+1))-1)) -1 ) << offset_len;
	} else {
		off = ranked_lz_offset[i] - 3 ;
		up_bitCount += prefix_len;
	}
	return (off);
}


uint32_t pucrunch_256_decode (uint8_t *outbuf,uint16_t *out_len,
                             const uint8_t *inbuf, uint16_t in_len) {

// start_clock();
	up_SetInput (inbuf);


	uint8_t alphabet[256];
	for (uint16_t i=0; i<256;i++) alphabet[i]=i;


	// decompressed length
	// datatype limits block size, needs to be signed for loop
	int32_t len = ((up_GetValue()-1)<<8) | up_GetBits(8);
	// number of ESCape bits and first ESCape sequence
	uint8_t lz_length_max_gamma = up_GetBits (4);
	uint8_t no_esc = up_GetBits(4);
	uint8_t current_esc8 = up_GetBits(no_esc);
	// RLE table count and table w/ up to 15 entries
	uint8_t rank_len = up_GetValueMaxGamma(6)-1;
	uint8_t rle_rank[MAX_RLE_RANKS];
	for (uint8_t i=0; i < rank_len; i++)
		rle_rank[i] = up_GetBits(8);
	// number of '1' bits to identify a following character (as opposed to rank pointer)
	uint8_t no_rle_char_esc_bits = (rank_len==0)?0 : int_log2(rank_len)+1;
	// number of bits for LZ2 offset
	uint8_t lz2_bits = up_GetBits(4);
	// number of LSBs for LZ offset
	uint32_t lz_lsb = up_GetBits(4);
	// LZ offset ranking
	uint32_t lz_rank_len_bits = up_GetValue()-1;
	uint32_t no_rank = (1 << lz_rank_len_bits)-1;
	uint32_t lz_rank_len = 0;
	uint32_t lz_rank_prefix_len = 0;
	uint32_t lz_rank_prefix;
	uint16_t lz_rank[127];
	if (lz_rank_len_bits > 0) {
		lz_rank_len = (1 << lz_rank_len_bits)-1;
		lz_rank_prefix_len = up_GetBits (4);
		lz_rank_prefix = up_GetBits (lz_rank_prefix_len);
		for (uint8_t i=0; i< lz_rank_len; i++) {
			lz_rank[i] = up_GetBits(lz_lsb);
			lz_rank[i] |= ((up_GetValue ()-1) << lz_lsb);
			lz_rank[i] += 3;
		}
	} else {

	}
	// move to front
	uint8_t mtf_second_line = up_GetValue()-1;
	// recursion level of LITerals' implicit huffman
	uint8_t lit_opt_huff_level = up_GetBits(3);
	// recursion level of LZ2 offsets' implicit huffman
	uint8_t lz2_opt_huff_level = up_GetBits(4);

//fprintf(stderr, "-----== DECODING STATS ==-----\n#ESC: 			%6u\n#LZ OFFSET LSB:		%6u\n#LZ2 OFFSET BITS	%6u\n#RANKED RLEs		%6u\nIN incl. header:	%6u\nOUT:			%6u\n------------------------------\n",
//	no_esc, lz_lsb, lz2_bits, rank_len, in_len, len);


	uint8_t rle1_c = 1<<no_rle_char_esc_bits;
	uint8_t rle8_c = 8-no_rle_char_esc_bits;
	uint8_t esc8_e = 8 - no_esc;
	uint8_t lzl_p  = lz_lsb - lz_rank_prefix_len;

	uint8_t current_esc = current_esc8 << (esc8_e);
	uint8_t esc_mask = 0xFF00 >> no_esc;

	for (uint16_t cnt=0;cnt<len;) {
		uint8_t first_fetch = up_GetBits(no_esc);
		if (likely(first_fetch == current_esc8)) {
//			uint16_t fetch = up_GetValueMaxGamma(lz_length_max_gamma);
			uint16_t fetch = up_GetValue();
			if (fetch == 1) { // '1' in E. Gamma means the first bit is '1'
				// investigate the following bit(s)
				switch (up_PeekBits(2)) {
					case 0:
						// RLE
						up_ConsumeBits (2);
						uint16_t count = up_GetValue()+1;
						uint8_t character = up_GetValueMaxGamma (no_rle_char_esc_bits);
						if (likely( character < (rle1_c)) ) {
							// ranked character
							character = rle_rank[character-1];
						} else {
							// not ranked, get rid of the preceding, and already read sequence of '0....01'
// only for character of type 16 bit:			character = (character-(rle1_c)) << (rle8_c);
							character <<= (rle8_c);
							// and get the remaining bits
							character |= up_GetBits(rle8_c);
						}
						// broadcast character
						uint64_t c = (character << 8) | character;
						c |= (c<<16);
						c |= (c<<32);
						uint64_t *d; d = &outbuf[cnt];
//						for (; count!=0; count--) outbuf[cnt++] = character;
						for (uint16_t i = 0; i < count; i+=sizeof(*d)) { *d=c; d++;}
						cnt += count;
						break;
					case 1:
						up_ConsumeBits(2);
						// ESCaped (maybe mtf-encoded) LITeral
						// get the new ESCape code
						current_esc8 = up_GetBits(no_esc);
						current_esc = current_esc8 << esc8_e;
						// fetch the rest of the mtf-encoded LITeral
						first_fetch <<= esc8_e;
					 	first_fetch |= up_get_huffed_value(esc8_e,lit_opt_huff_level);
						outbuf[cnt++] = move_to_front_decode (first_fetch, alphabet, mtf_second_line);
						break;
					default:
						// LZ - to be specific, this is LZ2
						up_ConsumeBits(1);
//						uint16_t offset = (lz2_bits<=2)?up_GetBits(lz2_bits):up_get_huffed_value (lz2_bits,lz2_opt_huff_level);
						uint16_t offset = up_get_huffed_value (lz2_bits,lz2_opt_huff_level);
						memcpy (&outbuf[cnt], &outbuf[cnt-2-offset], 2); cnt +=2;
						break;
				}
			} else {
				// common LZ, the fetch already contains the LZ count
				fetch++;
				uint16_t offset = up_get_lz_offset (lz_rank_prefix_len, lz_rank_prefix, lz_rank_len_bits, lz_lsb, lz_rank);
				// !!! requires 64K+sizeof(*dst) sized-buffer
				// !!! other platforms might run better on 32 bit
				uint64_t *src; src = &outbuf[cnt-3-offset];
				uint64_t *dst; dst = &outbuf[cnt];
				for (uint16_t i=0; i < fetch; i+=sizeof(*dst)) {
					*dst = *src; dst++; src++;
				}
				// memcpy (&outbuf[cnt], &outbuf[cnt-3-offset],fetch);
				cnt += fetch;
			}
		} else {
			// unESCaped LITeral, fetch the the rest
			first_fetch <<= esc8_e;
			first_fetch |= up_get_huffed_value(esc8_e,lit_opt_huff_level);
			outbuf[cnt++] = move_to_front_decode (first_fetch, alphabet, mtf_second_line);
		}
	}
	*out_len = len;
}



// ---=== MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN ===---

int main(int argc, char* argv[]) {

        pucr_init();
//    if(argc != 2)
//    {
//        printf( "usage:  %s [-f<{0,1,2,3}>]\ndata via stdin and stdout", argv[0] );
//        return 0;
//    }
//
        char* buffer;
        buffer=(char *)malloc(65536+1);
        size_t fileLen = fread(buffer, sizeof(char),65536,stdin);

        if (!buffer)
        {
                fprintf(stderr, "memory error!");
                return (-1);
        }

        fread(buffer, sizeof(char), fileLen, stdin);

        char output[65535];
        uint16_t outlen = 65535;

        pucrunch_256_encode (output, &outlen, buffer, fileLen, FAST_LANE);
        free(buffer);
        fprintf (stderr, "pucr: %u --> %u bytes\n",fileLen, outlen);
//      fwrite (out_buffer, sizeof(char), outlen, stdout);

print_clock(); start_clock();
for (uint16_t i=0; i < DECODE_ITERATIONS; i++) {
// uint16_t outlen_old = outlen;
        pucrunch_256_decode (output, &outlen, out_buffer, outlen);
// outlen=outlen_old;
}
double t = print_clock();
fprintf (stderr, "ITERATIONS: %u\n", DECODE_ITERATIONS);
fprintf (stderr, "--> %.1f MB/s\n", fileLen * DECODE_ITERATIONS / t / 1024 / 1024);

        fwrite (output, sizeof(char), outlen, stdout);

        return(0);
}
