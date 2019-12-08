#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>


#include <time.h>
static clock_t start, end;
static double cpu_time_used;

#define NO_MAX_GAMMA 0xFF

#define MAX_RLE_RANKS 31

#define LZ	0x0000
#define LIT	0x0002
#define RLE	0x0003
#define MARKED	0x0080


#define LZ_LEN_MAX_GAMMA 15
#define MAX_ESC 9
#define MAX_MTF 254

struct pucr_node {
	uint16_t	rle_count;
	uint16_t	lz_count;
	uint16_t	lz_max_offset;
	uint16_t	lz_cur_offset;
	uint32_t	lz_cost;
	uint16_t	seen_before;

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

static void start_clock(void) {
	start = clock();
}


static void print_clock(void) {
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	fprintf (stderr,"CPU TIME USED: %f\n", cpu_time_used);
}


static uint8_t int_log2_lut [65536];
static int8_t lz_offset_encoded_length[16][16];


static uint32_t pucr_init (void) {

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


static uint8_t int_log2 (uint16_t x) { // uint16_t x
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

/*struct alphabet_node {
	uint8_t character;
	uint8_t position;
};*/


static uint8_t move_to_front_encode (uint8_t in_char,
                                     uint8_t alphabet[],
				     uint8_t second_line) {

        // search letter in alphabet
        uint16_t j= 0;
        for (; j < 256; j++) if (alphabet[j]==in_char) break;
	uint8_t new = j>second_line?second_line:0;
        // move letter to front of alphabet
        for (uint16_t k=j; k>new; k--) alphabet[k]=alphabet[k-1];
        alphabet[new]=in_char;
        return(j);
}


static uint8_t move_to_front_decode (uint8_t in_char,
                                     uint8_t alphabet[],
				     uint8_t second_line) {

        uint8_t ret = alphabet[in_char];

        // move letter to front of alphabet
	uint8_t new = in_char>second_line?second_line:0;
        for (uint16_t k=in_char; k>new; k--) alphabet[k]=alphabet[k-1];
        alphabet[new]=ret;
        return(ret);
}


// ---=== LIT LIT LIT LIT LIT LIT LIT LIT LIT LIT LIT LIT LIT LIT LIT ===---


static uint32_t update_lit_occ_from_graph (struct pucr_node graph[], uint16_t in_len,
		   uint16_t lit_occ[]) {

	uint16_t i;
	for (i=0; i<256; i++) lit_occ[i] = 0;

	for (i=0; i < in_len; i=graph[i].next_node)
		if (graph[i].way_to_go == LIT)
			lit_occ[graph[i].cur_lit]++;

	return (0);
}


static int32_t huffiness (uint16_t occ[], uint16_t first, uint16_t last, uint8_t *opt_level) {

	// how 'huffy' are the values? does it make sense to encode them with an implicit huffman tree?
	uint16_t i;
	int32_t sum = 0;

	uint8_t huff_level = 0;
	int32_t max_sum = sum;
	uint16_t q;
	uint16_t j;
	for (i=last-first+1; i>=4; i=q) {
		q = i >> 2;
		for (j=first;j<(first+q);j++) sum += occ[j]; // one bit less
		for (j=first+(q<<1);j<first+i;j++) sum -= occ[j]; // one bit more
		if (sum > max_sum) {
			max_sum = sum;
			huff_level++;
		} else {
			break; // no need to go any further
		}
	}

	*opt_level = huff_level;
// fprintf (stderr, "FINAL BIT SUM: %i   OPT LEVEL: %u\n",max_sum, huff_level);
	return (max_sum);
}

// not used yet
static uint32_t estimate_lit_and_generate_occ_from_inbuf (struct pucr_node graph[], const uint8_t *inbuf, uint16_t in_len,
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
//				fprintf (stderr,"!!! RLE %c at %u !!!\n", inbuf[i+1],i);
			rle_count = 1;
			character = inbuf[i];
		}
		graph[i].rle_count = rle_count;

// fprintf (stderr, "Pos. %5u  LZ %5u  %c\n",i,rle_count,character); // !!!
	}
	// first position
	if (rle_count > 1) rle_occ[inbuf[i+1]]++;

// fprintf (stderr, "!!! RLE %c at %u !!!\n", inbuf[i+1],i+1); }
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
// fprintf (stderr, "Rnk. %2u Occ. %5ux [%3u]\n",i, rle_occ[max_index], max_index ); // !!!
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
// fprintf(stderr,"SUM RLE%u: %u\n",ranked,sum);
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
//				fprintf (stderr,"!!! RLE %c at %u !!!\n", inbuf[i+1],i);
			rle_count = 1;
			character = inbuf[i];
		}
		graph[i].rle_count = rle_count;

// fprintf (stderr, "Pos. %5u  RLE %5u  %c\n",i,rle_count,character); // !!!
	}
	// first position
	if (rle_count > 1) rle_occ[inbuf[i+1]]++;

// fprintf (stderr, "!!! RLE %c at %u !!!\n", inbuf[i+1],i+1); }
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
// fprintf (stderr, "   * closer match at offset %u (offset: %u instead of %u)\n",pattern_start,graph[pattern_start].lz_cur_offset,graph[pattern_start].lz_max_offset);
	                break;
                }
        }
}


static uint32_t prepare_last_seen_from_graph (struct pucr_node graph[], const uint8_t *inbuf, uint16_t in_len) {

	uint16_t last_occ [65536];
	memset (last_occ, 0xFF, 0x20000);

	for (uint16_t i=0; i < in_len;i++ ) {
// if (inbuf[i] == inbuf[i+1]) continue;
		uint16_t idx = (inbuf[i] << 8) | inbuf[i+1];
		graph[i].seen_before = last_occ[idx];
		last_occ[idx] = i;

	}
}


static uint32_t find_matches (struct pucr_node graph[], const uint8_t *inbuf, uint16_t in_len,
				uint16_t first_pos, uint16_t last_pos,
				uint8_t recursive, uint8_t fast_lane, uint8_t lz_length_max_gamma) {

	// a block should not be much shorter than 256...512 bytes
	if ( (recursive > 0) && ((last_pos - first_pos) > 512) ) {
		uint16_t sep = (last_pos - first_pos) * 43 / 64;
		find_matches (graph, inbuf, in_len, first_pos, first_pos+sep, recursive-1, fast_lane, lz_length_max_gamma);
		find_matches (graph, inbuf, in_len, first_pos+sep+1, last_pos, recursive-1, fast_lane, lz_length_max_gamma);
	} else {
		uint16_t i,j,k, off;
		for (i=first_pos; i <= last_pos;i+=fast_lane?graph[i].rle_count:1) {
			graph[i].lz_max_offset = 0xFFFF; // !!! neccessary for msb-lsb-offset calculation of number of bits
			graph[i].lz_count =  0x0000; // = no  match yet
graph[i].lz_cost = 0;
			// we will be able to skip at least 2 characters when comparing
			// because we know we already have a 2-byte match.
			// even more if rle is longer in that position
			off = (graph[i].rle_count -1);
			off = off<=2?2:off;
// off = 0;
			for (j=graph[i].seen_before; j != 0xFFFF; j=graph[j].seen_before) { // check at offset J for a match ...
				if (graph[i].rle_count != graph[j].rle_count) continue;
				k = j + off;
				// if a match crosses last_pos, "< in_len" would be advantageous to "<=last_pos" ... !!! may be consider for fast_lane
				for (k=j+off; (k < i) && ((i+(k-j))<in_len) && ( (k-j) < (3<<(lz_length_max_gamma-1)) ); k++) { // ... of length (K-J)
					if (inbuf[k] != inbuf[i+(k-j)]) break;
				}
				if (k > i) k = j ; // no match
				else  {
					k -= j; // match length
					// treatment of a match found
					// naive weighting function: the longer the better // replace with weighting/boundary conditions later on !!!
// experimental: try to calculate cost.
// WHY does it result in a compression gain?!
uint32_t gain = (k==2)?16-(3+8):
	               k*8-(1 + int_log2(k-1)*2+1 + (i-j-k)<1024?11:10+int_log2((i-j-k)>>10)*2+1);


//					if (k > graph[i].lz_count) {
					if (gain > graph[i].lz_cost) {
							graph[i].lz_max_offset = i-j; // store the 'real' offset
							graph[i].lz_count=k; // ; -(off2-i); // k
graph[i].lz_cost = gain;
					}
				}
			}
		}
fprintf (stderr, "--- FINDING MATCHES %5u ... %5u: ",first_pos,last_pos);
print_clock(); start_clock();
	}
}



int compare( const void* a, const void* b)
{
return ( *(int*)a - *(int*)b );
//     return (a < b) - (a > b);
}




static uint32_t find_optimal_lz_offset_no_lsb (struct pucr_node graph[], uint16_t in_len,
					       uint8_t no_esc,
					       uint8_t *no_lz_offset_lsb, uint8_t *no_lz2_offset_bits,
					       uint8_t *lz2_opt_huff_level) {

	uint16_t lz_occ[16] = {0};
	uint16_t lz2_occ[16] = {0};

uint32_t lz_offset[65536]= {0};

uint16_t lz2_full_occ[65536] = {0};

	// after all matches are finally determined, calculate the optimum of lsb-msb coding ratio
	// first step: count them, especially the bits their offsets need
	for (uint16_t i=0; i != in_len; i=graph[i].next_node) {
		if (graph[i].way_to_go == LZ) {
			// lz_length is not necessarily the real way to go,
			// but next_node is
lz_offset[graph[i].lz_cur_offset]++;
			if ((graph[i].next_node - i) >= 3) {
				// for the 3-byte and longer matches: calculate the case with "p" pure LSB and the rest E.Gamma-coded
			        for (int p=0; p <16; p++) // 0 = 1-bit offset, 1 = 2-bit offset, ...
					lz_occ[p] += (p+1) + int_log2( ((graph[i].lz_cur_offset-(graph[i].next_node-i))>>(p+1))+1 )*2+1;
			} else {
				lz2_full_occ[graph[i].lz_cur_offset-2]++;
				int p = int_log2 (graph[i].lz_cur_offset - 2);
				p = p<=15?p:15;
				lz2_occ[p]++;
			}
		}
	}

// qsort( lz_offset, 65536, sizeof(uint32_t), compare );
// for (uint16_t g=0; g < 65535;g++) fprintf (stderr, "LZ OFFSET Rnk %u. %u x\n",g,lz_offset[g]);
// for (uint16_t g=0; g < 16;g++) fprintf (stderr, "LZ OFFSET Rnk %u. %u x\n",g,lz_offset[g]);

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
// fprintf (stderr, "LZ2 OFFSET %u BITS SAVED %d   LOST %d   RESULT: %d [maxp: %u]\n",p+1, gained, lost, gained-lost, max_p);
	}
// fprintf (stderr, "ESC: %2u LZ2 OFFSET %u BITS\n",no_esc, min_lz2_bits);
	*no_lz2_offset_bits = min_lz2_bits;
uint32_t lz2_huff_savings = huffiness (lz2_full_occ,0,(1<<(min_lz2_bits))-1, lz2_opt_huff_level);
// fprintf (stderr, "LZ2 SAVINGS: %u     OPT HUFF: %u   (LZ MIN: %u)\n",lz2_huff_savings,*lz2_opt_huff_level,min_lz2_bits);

	return (0);
}


// ---=== ESC ESC ESC ESC ESC ESC ESC ESC ESC ESC ESC ESC ESC ESC ESC ===---

static uint8_t get_huffed_prefix (uint8_t value, uint8_t bits, uint8_t prefix_len) {

	uint16_t mask = 1<<(bits-1);
	mask |= mask >> 1;

	if ((value & mask) == 0)
		return (value >> (bits-1-prefix_len));
	else if ((value & mask) == (1<<bits-2)) {
//		put_bit (1);
//		put_n_bits (value, bits-1);
		return ((1 << (prefix_len-1))  | (value >> (bits-prefix_len)));
	} else {
//		put_bit (1);
//		put_bit (0);
//		put_n_bits (value, bits-1);
		return ((1 << (prefix_len-1))  | (value >> (bits-prefix_len+1)));

	}
}




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
            graph[i].cur_lit = mtf?move_to_front_encode (inbuf[i], alphabet, mtf_second_line):inbuf[i];
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
		uint32_t rle_cost, lz_cost, lit_cost;
		uint32_t min = 0xFFFFFFFF;
		uint16_t rle_way, lz_way;

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
				  int_log2(j-1)*2+1 + lz_lsb + int_log2(((graph[i].lz_max_offset-j)>>lz_lsb)+1)*2+1;
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
						  int_log2(used_match_length-1)*2+1 + lz_lsb + int_log2(((graph[i].lz_cur_offset-used_match_length)>>lz_lsb)+1)*2+1;
					lz_cost += graph[i+used_match_length].bits_to_end+no_esc;
				}
			}
		} else {
			lz_cost = min;

		}
		/* determine LIT cost (always possible), naively add 1 for
		    ocasionally occuring LIT code "010" and new ESC bits.

		    This will not lead to the correct number of bits but takes
		    them into account for decision. Their real number is
		    counted in the optimzize_escape step */
		lit_cost = 8+graph[i+1].bits_to_end ; // + 1 ??? !!! ;

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
}


// ---=== OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT ===---


static uint16_t outPointer = 0;
static uint8_t bitMask = 0x80;
static uint32_t bitCount = 0;
#define OUT_SIZE 65536
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

    while (value>1) {
        bits = (bits<<1) | (value & 1); /* is reversed compared to value */
        value >>= 1;
        count++;
        put_bit (1);
    }
    if (count<maxGamma)
        put_bit (0);
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
			uint8_t lz_lsb, uint8_t lz2_bits, uint8_t lz_length_max_gamma, uint8_t lz2_opt_huff_level) {

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
	// move to front
	if (mtf) put_bit (1); else put_bit (0);
	if (mtf) put_n_bits (mtf_second_line,8);

	// re-think sharply: the number of bits to be encoded is limited by the number of esc-bits (and thus remaining lit-bits) !!!
	if (mtf) put_n_bits (lit_opt_huff_level, 3); // could reach up to 4 on case of 0 esc bits

	// re-think sharply: the number of bits to be encoded is limited by the number of lz2 offset-bits !!!
	put_n_bits (lz2_opt_huff_level,4);
	// derive further parameters
	uint8_t no_rle_char_esc_bits = (rle_rank_len==0)?0 : int_log2(rle_rank_len)+1;
// fprintf (stderr, "RLE CHAR ESC CODE: no of 1's: %u\n", no_rle_char_esc_bits);

int16_t old_lz_offset[4] = {0,0,0,0};
uint16_t esc_seq_count = 0;
uint16_t unesc_seq_count = 0;
uint16_t no_unesc_seq = 0;
uint16_t sum_unesc_count = 0;
uint16_t lz_off[256] = {0};
uint16_t lz_off_hist_match = 0;
uint16_t rle_len[256] = {0};
uint16_t missed_rle2  = 0;
uint16_t rle2_ranked = 0;
uint16_t rle2_unranked = 0;
uint8_t rle_char_hist = 0;
uint16_t rle_char_hist_match = 0;
	// output data
	for (uint16_t i=0; i!=in_len; i=graph[i].next_node) {
// int old = bitCount;
		switch (graph[i].way_to_go) {
			case RLE:
fprintf(stderr,"");
uint16_t rle = graph[i].next_node -i;
if (inbuf[i]==rle_char_hist) {
rle_char_hist_match++;
}

rle_char_hist=inbuf[i];


rle = rle>63?63:rle;
rle_len[rle]++;

				// ESC ++ 011
				put_n_bits (current_esc_8, no_esc);
				put_bit (0); put_bit (1); put_bit (1);
                                // output RLE count (real count minus 1, i.e. 1 instead of 2 etc.) E. Gamma coded
                                put_value (graph[i].next_node - i - 1, NO_MAX_GAMMA);
				// check if it is a ranked character
				uint8_t k;
				for (k=0; k < rle_rank_len; k++) if ( inbuf[i] == rle_rank[k] ) break;
				if (k != rle_rank_len) {
					// is ranked: E. Gamma coded pointer to rank list
					put_value (k+1, NO_MAX_GAMMA);
if (rle==2) rle2_ranked++;
//					fprintf (stderr, " {%u} ", no_esc + 3 + int_log2(k+1)*2+1 + int_log2(graph[i].next_node - i - 1)*2+1);
				} else {
					// not ranked: determine '11..11'-code length (already determined outside loop) followed by character
					put_n_bits (0xFF, no_rle_char_esc_bits);
					put_n_bits (inbuf[i], 8);
if (rle==2) rle2_unranked++;
//					fprintf (stderr, " {%u} ", no_esc + 3 + no_rle_char_esc_bits + 8 + int_log2(graph[i].next_node - i - 1)*2+1);
				}
//				fprintf (stderr, "%u. RLE(%c,%u) --> %u [%u vs %u]\n",i,inbuf[i],graph[i].next_node-i,graph[i].next_node,graph[i].bits_to_end-graph[graph[i].next_node].bits_to_end, bitCount-old);

esc_seq_count++;
if (unesc_seq_count > 8) { fprintf (stderr, "UNESC SEQ: %u\n", unesc_seq_count);
sum_unesc_count += unesc_seq_count; no_unesc_seq++; }
unesc_seq_count=0;

				break;
			case LZ:
/* if ( (graph[i].next_node-i) == 2 )	fprintf (stderr, "%u. LZ(%u,%u) --> %u [%u]\n",i,
										graph[i].next_node-i,
										graph[i].lz_cur_offset,
										graph[i].next_node,
										-graph[graph[i].next_node].bits_to_end + graph[i].bits_to_end);
*/				put_n_bits (current_esc_8, no_esc);

uint16_t lz= graph[i].lz_cur_offset-(graph[i].next_node -i);
// lz = lz>255?255:lz;
// if ((graph[i].next_node-i)!=2) lz_off[lz]++;

				if ( (graph[i].next_node-i) == 2 ) {
					// LZ2: ESC + 000 + LZ2 OFFSET (lz2_bits)
	                                put_bit (0); put_bit (0); // put_bit (0);
if (lz2_bits <= 2)
					put_n_bits (graph[i].lz_cur_offset-2, lz2_bits);
else
put_huffed_value (graph[i].lz_cur_offset-2, lz2_bits,lz2_opt_huff_level);
//					fprintf (stderr, " {%u} ", no_esc + 2 + lz2_bits);
				} else {
uint16_t n=4; // history length, we use 1 for now !!!
uint16_t j=0;
for (;j<n;j++)
if (lz == old_lz_offset[j]) {
	lz_off_hist_match++;
	break;
}

if (j==n) {
	old_lz_offset[3] = old_lz_offset[2];
	old_lz_offset[2] = old_lz_offset[1];
	old_lz_offset[1] = old_lz_offset[0];
	old_lz_offset[0] = lz;
} else {
	for (uint16_t k=j;k>0; k--)
		old_lz_offset[k] = old_lz_offset[k-1];
	old_lz_offset[0] = lz;
}
					// LZ: ESC + E. Gamma coded actually used lz length (next_node - i - 1), not neccessariliy lz_max_len ...
//					put_value (graph[i].next_node - i - 1 , lz_length_max_gamma);
j = graph[i].next_node-i-1;
put_value (j , lz_length_max_gamma);
					// ... and also E. Gamma coded LZ offset minus the actually used length of lz match(=next_node-i)
					uint16_t off = (graph[i].lz_cur_offset-(graph[i].next_node-i));
					// MSBs (+1 for E. Gamma code)
					put_value ((off >> lz_lsb)+1, NO_MAX_GAMMA);
					// LSBs
					put_n_bits (off, lz_lsb);
// put_huffed_value (off, lz_lsb,0);
//					fprintf (stderr, " {%u} ", no_esc + int_log2(graph[i].next_node-i-1)*2+1 + int_log2((off>>lz_lsb)+1)*2+1 + lz_lsb);
				}
/* fprintf (stderr, "%u. LZ(%u,%u) --> %u [%u vs bitCount %u]\n",i,
										graph[i].next_node-i,
										graph[i].lz_max_offset,
										graph[i].next_node,
										-graph[graph[i].next_node].bits_to_end + graph[i].bits_to_end, bitCount-old);
*/
esc_seq_count++;
if (unesc_seq_count > 8) { fprintf (stderr, "UNESC SEQ: %u\n", unesc_seq_count);
sum_unesc_count += unesc_seq_count; no_unesc_seq++; }
unesc_seq_count=0;
				break;

			case LIT:
if (inbuf[i] == inbuf[graph[i].next_node] && graph[graph[i].next_node].way_to_go == LIT) 
missed_rle2++;
//				fprintf (stderr, "%u. LIT(%c) --> %u [%u]",i,inbuf[i],graph[i].next_node, -graph[graph[i].next_node].bits_to_end + graph[i].bits_to_end+(graph_current_esc == graph[i].new_esc?0:no_esc+3));

				if (graph[i].new_esc != graph_current_esc) {
//					fprintf (stderr, "  NEW possible ESC: %u", graph[i].new_esc);
					graph_current_esc = graph[i].new_esc;
				}
				if ( (graph[i].cur_lit & current_esc_mask) != current_esc) {

unesc_seq_count++;
if (esc_seq_count > 8) fprintf (stderr, "ESC SEQ: %u\n", esc_seq_count);
esc_seq_count=0;

					// output the mtf-encoded literal
// normal:					put_n_bits (graph[i].cur_lit, 8);
//put_huffed_value(graph[i].cur_lit,8);

if ( mtf == 0 || no_esc > 6 || graph[i].cur_lit > 0xFF>>no_esc) {
put_n_bits (graph[i].cur_lit, 8);
} else {
put_n_bits (graph[i].cur_lit >> (8-no_esc), no_esc);
put_huffed_value(graph[i].cur_lit, (8-no_esc),lit_opt_huff_level);
} 					// select huffman tree if available
					// eventually, use codebook for output
/*					if ( head_array[(graph[i].cur_lit) >>(8-no_esc)] == 0xFFFF)
						put_n_bits (graph[i].cur_lit, 8- no_esc);
					else {
						for (uint8_t j=0; j < table[graph[i].cur_lit].code.code_len; j++) {
							if ( table[graph[i].cur_lit].code.code[j >> 3] & (0x80 >> (j & 7)) )
								put_bit(1);
							else
								put_bit(0);
						}

					}
*/
// fprintf (stderr, " {%u} ",8);
				} else {
esc_seq_count++;
if (unesc_seq_count > 8) { fprintf (stderr, "UNESC SEQ: %u\n", unesc_seq_count);
sum_unesc_count += unesc_seq_count; no_unesc_seq++; }
unesc_seq_count=0;


					// ESC = first bits of LITeral
					put_n_bits (current_esc_8, no_esc);
					// 010
	                                put_bit (0); put_bit (1); put_bit (0);
					// new ESC code
					current_esc = graph[i].new_esc;
					current_esc_8 = current_esc >> (8 - no_esc);
					put_n_bits (current_esc_8, no_esc);
//					fprintf (stderr, "NEW ESCAPE: %u ", current_esc);
					// rest of mtf-encoded LITeral:
					// select huffman tree if available
					// eventually, use codebook for output

// normal:					put_n_bits (graph[i].cur_lit, 8-no_esc);
if ( mtf == 0 || no_esc > 6 || graph[i].cur_lit > 0xFF>>no_esc)
put_n_bits (graph[i].cur_lit, 8-no_esc);
else {
put_huffed_value(graph[i].cur_lit, (8-no_esc),lit_opt_huff_level);
}
/*					if ( head_array[(graph[i].cur_lit) >>(8-no_esc)] == 0xFFFF){
						put_n_bits (graph[i].cur_lit, 8- no_esc);
					} else {
	                                        for (uint8_t j=0; j < table[graph[i].cur_lit].code.code_len; j++) {
	                                                if ( table[graph[i].cur_lit].code.code[j >> 3] & (0x80 >> (j & 7)) )
	                                                        put_bit(1); else put_bit(0);
	                                        }
					}
*/
// fprintf (stderr, " {%u} ",8+3+no_esc);
				}
				break;
			default:
				// !!! shall not happen
				break;
		}
// fprintf (stderr,"\n");
	}


//	printf ("OUT POINTER: %u  BITMASK: %u\n", outPointer, bitMask);
	flush_bits();
	*out_len = (bitCount+7) >> 3;
// fprintf (stderr, "TOTAL SIZE: %u BITS ~~ %u BYTES\n", bitCount, (bitCount+7) >> 3);
//	printf ("OUT POINTER: %u\n", outPointer);
// for(uint16_t h=0; h< 256;h++) fprintf(stderr, "LZ OFF%3u: %5u x    RLE%3u: %5u\n",h,lz_off[h],h,rle_len[h]);
fprintf (stderr, "SUM UNESC SEQ COUNT: %u (%u x)\n", sum_unesc_count, no_unesc_seq);
fprintf (stderr, "missed RLE2 opportunities: %u\n", missed_rle2);
fprintf (stderr, "unranked RLE2s: %u\n", rle2_unranked);
fprintf (stderr, "ranked RLE2s: %u\n", rle2_ranked);
fprintf (stderr, "RLE CHAR HISTORY MATCH: %u x\n", rle_char_hist_match);
fprintf (stderr, "LZ OFFSET HISTORY MATCH: %u x\n", lz_off_hist_match);

}


void print_graph_stats (const struct pucr_node graph[], uint16_t in_len, uint8_t no_esc, uint8_t start_esc) {

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
	uint16_t lz_bytes = 0;

	uint8_t current_esc = start_esc;
	uint8_t graph_current_esc = start_esc;

	uint8_t current_esc_mask = 0xFF00 >> no_esc;

	uint8_t min_esc_seq_len = 16 >> no_esc;

	for (i=0; i != in_len; i=graph[i].next_node) {
		uint8_t esc = 1;
		uint8_t lit = 0;
		if ( graph[i].way_to_go == RLE) {
			rle_cnt++;
			rle_bytes += graph[i].next_node - i;
		}
		if ( graph[i].way_to_go == LZ) {
			lz_cnt++;
			lz_bytes += graph[i].next_node - i;
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
		} else { // any non-LIT

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
#LZ CNT:		%u\n\
LZed BYTES:             %u\n\
------------------------------\n",
rle_cnt, rle_bytes, lit_cnt,lit_esc_cnt,lit_unesc_cnt, unesc_seq_cnt,unesc_seq_lit, min_esc_seq_len,esc_seq_cnt,esc_seq_x, lz_cnt,lz_bytes);

for (i=0; i< 256; i++) {
// fprintf (stderr, "LIT %3u: %5u x\n",i, cur_lit_occ[i]);
}

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
	uint8_t lit_rank[256];
	uint8_t lit_rank_len;

	count_rle_and_generate_occ_from_inbuf (graph, inbuf, in_len, rle_occ);
	generate_rle_ranks (rle_occ, rle_rank, &rle_rank_len);
fprintf (stderr, "--- COUNTING RLE: ");
print_clock(); start_clock();
	prepare_last_seen_from_graph (graph, inbuf, in_len);
fprintf (stderr, "--- PREPARE OCCURENCES: ");
print_clock(); start_clock();

uint8_t lz_length_max_gamma = LZ_LEN_MAX_GAMMA;

	find_matches (graph, inbuf, in_len, 0, in_len-1, 2, fast_lane, lz_length_max_gamma);

// fprintf (stderr, "--- FINDING MATCHES TOTAL: ");
// print_clock(); start_clock();
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

			find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rle_rank_len, no_esc, lz_lsb, lz2_bits, fast_lane);
// fprintf (stderr, "AFTER FIRST OPTIMIZE: %u bits\n",graph[0].bits_to_end);
			find_optimal_lz_offset_no_lsb (graph, in_len, no_esc, &lz_lsb, &lz2_bits, &lz2_opt_huff_level);

// too early determination seems to lead to unreliable results...
// thus, a generic value is chosen here, the real one chosen
// futher down the loop
if (fast_lane < 2) lz_lsb = 10;

// fprintf (stderr, "NO ESC: %u --- LZ2: %u --- LZLSB: %u \n",no_esc,lz2_bits,lz_lsb);
			int escaped = optimize_esc (graph, inbuf, in_len, no_esc, &start_esc, mtf, 0);

if (fast_lane < 2)	update_rle_occ_from_graph (graph, inbuf, in_len, rle_occ);
if (fast_lane < 2)	generate_rle_ranks (rle_occ, rle_rank, &rle_rank_len);
if (fast_lane < 2)	find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rle_rank_len, no_esc, lz_lsb, lz2_bits, fast_lane);
// fprintf (stderr, "AFTER SECOND OPTIMIZE: %u bits\n",graph[0].bits_to_end);
if (fast_lane < 2)	find_optimal_lz_offset_no_lsb (graph, in_len, no_esc, &lz_lsb, &lz2_bits, &lz2_opt_huff_level);

	            	/* Compare value: bits lost for escaping -- bits lost for prefix */
			if ((graph[0].bits_to_end + (no_esc+3) * escaped) < min_no_bits) {
				min_no_esc = no_esc;
				min_no_bits = graph[0].bits_to_end + (no_esc+3)*escaped; // "- escaped" is for those counted 9-bit LIT in find_best_path
				// also store the coresponding lz_lsb value
				min_lz_lsb = lz_lsb;
				min_lz2_bits = lz2_bits;
				min_lz2_opt_huff_level = lz2_opt_huff_level;
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
	}
	// (re-)calculate best path with the values found or set above
	// alternatively, store the best graph found in the loop so far... (memory!!!)
	find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rle_rank_len, min_no_esc, min_lz_lsb, min_lz2_bits, fast_lane);
if (fast_lane < 2)		update_rle_occ_from_graph (graph, inbuf, in_len, rle_occ);
if (fast_lane < 2)		generate_rle_ranks (rle_occ, rle_rank, &rle_rank_len);
if (fast_lane < 2)		find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rle_rank_len, min_no_esc, min_lz_lsb, min_lz2_bits, fast_lane);
//if (fast_lane < 2)		update_rle_occ_from_graph (graph, inbuf, in_len, rle_occ);
//if (fast_lane < 2)		generate_rle_ranks (rle_occ, rle_rank, &rle_rank_len);


	// handling of remaining LITerals:
	int escaped;
	int32_t huffed = 0;

	uint8_t opt_level;

fprintf (stderr, "GENERATING FINAL GRAPH: ");
print_clock(); start_clock();

	// determine if we would be better off with or without mtf
	uint8_t min_mtf = 0;
if (fast_lane < 2) {
	uint32_t min_mtf_bits = 0xFFFFFFFF;
	for (mtf = 0; mtf < 2; mtf++) {
		escaped = optimize_esc (graph, inbuf, in_len, min_no_esc, &min_start_esc, mtf, 0);
		min_no_bits = graph[0].bits_to_end + (min_no_esc+3)*escaped; // "- escaped" is for those counted 9-bit LIT in find_best_path
huffed = 0;
if (!fast_lane)		update_lit_occ_from_graph (graph, in_len, lit_occ);
if (!fast_lane)		huffed = huffiness (lit_occ, 0, 0xFF>>min_no_esc, &opt_level);
// fprintf (stderr, "min_no_bits = %u\n",min_no_bits);
		if ( (min_no_bits-huffed) < min_mtf_bits ) {
			min_mtf = mtf;
			min_mtf_bits = (min_no_bits - huffed);
		}
	}
}
	// if mtf, will a higher mtf_second_line be beneficial?
	uint8_t min_mtf_second_line = 0;
	uint8_t  min_lit_opt_huff_level = 0;
	uint32_t min_mtf_second_line_bits = 0xFFFFFFFF;
	if (min_mtf) {
		for (uint8_t mtf_second_line=0; mtf_second_line < (fast_lane?1:MAX_MTF); mtf_second_line++) {
			escaped = optimize_esc (graph, inbuf, in_len, min_no_esc, &min_start_esc,min_mtf, mtf_second_line);
			min_no_bits = graph[0].bits_to_end + (min_no_esc+3)*escaped;
huffed=0;
update_lit_occ_from_graph (graph, in_len, lit_occ);
huffed = huffiness (lit_occ, 0, 0xFF >> min_no_esc, &opt_level);
// fprintf (stderr, "MTF CHECK %u BITS: %u\n",mtf_second_line, min_no_bits-huffed);

			if ( (min_no_bits-huffed) < min_mtf_second_line_bits ) {
// fprintf (stderr,"huffed = %u\n",huffed);
				min_mtf_second_line = mtf_second_line;
				min_mtf_second_line_bits = (min_no_bits - huffed);
				min_lit_opt_huff_level = opt_level;
			}
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
	output_path (graph, inbuf, in_len, out_len, rle_occ, rle_rank, rle_rank_len, min_start_esc, min_no_esc, min_mtf, min_mtf_second_line, min_lit_opt_huff_level, min_lz_lsb, min_lz2_bits, lz_length_max_gamma, min_lz2_opt_huff_level);

// fprintf (stderr, "START ESC: %u\n",min_start_esc);


fprintf(stderr, "-----== ENCODING STATS ==-----\n#ESC: 			%u\n#ESCAPED LIT:		%u\n#LZ OFFSET LSB:		%u\n#LZ2 OFFSET BITS	%u\n#RANKED RLEs		%u\nIN:			%u\nOUT: without header	%u\n     before huffman\n------------------------------\n",
	min_no_esc, escaped, min_lz_lsb, min_lz2_bits, rle_rank_len, in_len, (min_no_bits + 7)>>3);
print_graph_stats (graph, in_len,min_no_esc, min_start_esc);


fprintf (stderr, "--- OUTPUT FINAL PATH: ");
print_clock(); start_clock();
}

// ---=== DECODE DECODE DECODE DECODE DECODE ===---

static const unsigned char *up_Data;
static int up_Mask, up_Byte;
static uint32_t up_bitCount;

static void up_SetInput(const unsigned char *data) {
    up_Data = data;
    up_Mask = 0x80;
    up_Byte = 0;
    up_bitCount = 0;
}


static int up_GetBits(int bits) {
    int val = 0;

    while (bits--) {
	up_bitCount++;
        val <<= 1;
        if ((*up_Data & up_Mask))
           val |= 1;
        up_Mask >>= 1;
        if (!up_Mask) {
            up_Mask = 0x80;
            up_Data++;
            up_Byte++;
        }
    }
    return val;
}


static int up_GetValue(void) {
    int i = 0;

    while (1) {
        if (!up_GetBits(1))
            break;
        i++;
    }
    return (1<<i) | up_GetBits(i);
}


static int up_GetValueMaxGamma(int maxGamma) {
    int i = 0;

    while (i<maxGamma) {
        if (!up_GetBits(1))
            break;
        i++;
    }
    return (1<<i) | up_GetBits(i);
}


static int up_get_huffed_value (int bits, uint8_t rec) {

	if (rec == 0) return (up_GetBits(bits));

	if (!up_GetBits(1))
		if (rec == 1)
			return (up_GetBits(bits-2));
		else {
			return (up_get_huffed_value (bits-2,rec-1));
		}
	else if (up_GetBits(1))
		return (up_GetBits(bits-2) | (1<<(bits-2)));
	else
		return (up_GetBits(bits-1) | (1<<(bits-1)));
}


static uint32_t pucrunch_256_decode (uint8_t *outbuf,uint16_t *out_len,
                             const uint8_t *inbuf, uint16_t in_len) {

start_clock();
	up_SetInput (inbuf);

	uint8_t alphabet[256]; for (uint16_t i=0; i<256;i++) alphabet[i]=i;

	// decompressed length
	// datatype limits block size, needs to be signed for loop
	int32_t len = (up_GetValue()-1)<<8;
	len += up_GetBits(8);
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
	uint8_t lz_lsb = up_GetBits(4);
	// move to front
	uint8_t mtf = up_GetBits(1);
	uint8_t mtf_second_line = mtf?up_GetBits(8):0;
	// recursion level of LITerals' implicit huffman
        uint8_t lit_opt_huff_level = mtf?up_GetBits(3):0;
	// recursion level of LZ2 offsets' implicit huffman
	uint8_t lz2_opt_huff_level = up_GetBits(4);

fprintf(stderr, "-----== DECODING STATS ==-----\n#ESC: 			%6u\n#LZ OFFSET LSB:		%6u\n#LZ2 OFFSET BITS	%6u\n#RANKED RLEs		%6u\nIN incl. header:	%6u\nOUT:			%6u\n------------------------------\n",
	no_esc, lz_lsb, lz2_bits, rank_len, in_len, len);

	uint8_t current_esc = current_esc8 << (8 - no_esc);
	uint8_t esc_mask = 0xFF00 >> no_esc;
uint16_t rle_ranked_cnt = 0;

uint16_t c=0;
uint16_t old_lz_offset[4] = { 0 };

	for (uint16_t cnt=0;cnt<len;) {
		uint8_t first_fetch = up_GetBits(no_esc);
		if (first_fetch == current_esc8) {
			uint16_t fetch = up_GetValueMaxGamma(lz_length_max_gamma);
			if (fetch == 1) { // '1' in modified E. Gamma means the first bit is '0'
				// investigate the following bit(s)
				if ( up_GetBits(1) ) {
					if ( up_GetBits(1) ) {
						// RLE
						uint16_t count = up_GetValue()+1;
						uint16_t character = up_GetValueMaxGamma (no_rle_char_esc_bits);
						if ( character < (1<<no_rle_char_esc_bits) ) {
							// ranked character
							character = rle_rank[character-1];
rle_ranked_cnt++;
						} else {
							// not ranked, get rid of the preceding, and already read sequence of '1....1'
							character = (character-(1<<no_rle_char_esc_bits)) << (8-no_rle_char_esc_bits);
							// and get the remaining bits
							character |= up_GetBits(8-no_rle_char_esc_bits);
						}
//fprintf (stderr, "%u. RLE (%u,%c)\n", cnt, count,character);
						for (; count!=0; count--) outbuf[cnt++] = character;
					} else {
						// ESCaped (maybe mtf-encoded) LITeral
						// get the new ESCape code
						current_esc8 = up_GetBits(no_esc);
						current_esc = current_esc8 << (8 - no_esc);
						// fetch the rest of the mtf-encoded LITeral
// normal						first_fetch = (first_fetch<<(8-no_esc)) | up_GetBits(8-no_esc);
if (!first_fetch)					first_fetch = (first_fetch<<(8-no_esc)) | (mtf==1 && no_esc<=6?up_get_huffed_value(8-no_esc,lit_opt_huff_level):up_GetBits(8-no_esc));
else first_fetch = (first_fetch<<(8-no_esc)) | up_GetBits(8-no_esc);

						outbuf[cnt++] = mtf?move_to_front_decode (first_fetch, alphabet, mtf_second_line):first_fetch;
//fprintf (stderr, "%u. LIT (%c)  --ESCAPED--\n", cnt-1, first_fetch);
					}
				} else {
					// LZ - to be specific, this is LZ2
// normal					uint16_t offset = up_GetBits (lz2_bits);
					uint16_t offset = (lz2_bits<=2?up_GetBits(lz2_bits):up_get_huffed_value (lz2_bits,lz2_opt_huff_level));
// fprintf (stderr, "%u. LZ (2,%u)\n", cnt, offset+2);
					outbuf[cnt]=outbuf[cnt-2-offset]; cnt++;
					outbuf[cnt]=outbuf[cnt-2-offset]; cnt++;
				}
			} else {
				// common LZ, the fetch already contains the LZ count

uint16_t n=0;

if (fetch<(n+2)) {
	c++;
if (c<5)	fprintf (stderr, "%5u. len=%5u [j=%u]\n", c,old_lz_offset[fetch],fetch);
	uint16_t k=fetch-2;
	fetch = old_lz_offset[k];
	for (;k>0; k--)
		old_lz_offset[k] = old_lz_offset[k-1];
	old_lz_offset[0] = fetch;
} else {
	old_lz_offset[3] = old_lz_offset[2];
	old_lz_offset[2] = old_lz_offset[1];
	old_lz_offset[1] = old_lz_offset[0];
	fetch = fetch-n;
	old_lz_offset[0] = fetch;
}
				fetch++;
				// E. Gamma coded offset MSBs (-1)
				uint16_t offset = up_GetValue()-1;
				offset <<= lz_lsb;
				// LSBs
				offset |= up_GetBits (lz_lsb);
//fprintf (stderr, "%u. LZ (%u,%u)\n", cnt, fetch, offset+fetch);
				for (uint16_t i=0; i < fetch; i++) {
					outbuf[cnt]=outbuf[cnt-fetch-offset]; cnt++;
				}
			}
		} else {
			// unESCaped LITeral, fetch the the rest
// normal			first_fetch = (first_fetch<<(8-no_esc)) | up_GetBits(8-no_esc);
if (!first_fetch)	first_fetch = (first_fetch<<(8-no_esc)) | (mtf==1 && no_esc<=6?up_get_huffed_value(8-no_esc,lit_opt_huff_level):up_GetBits(8-no_esc));
else first_fetch = (first_fetch<<(8-no_esc)) | up_GetBits(8-no_esc);

			outbuf[cnt++] = mtf?move_to_front_decode (first_fetch, alphabet, mtf_second_line):first_fetch;
//fprintf (stderr, "%u. LIT (%c)\n", cnt-1, first_fetch);
		}
	}
print_clock(); start_clock();
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

	pucrunch_256_encode (output, &outlen, buffer, fileLen, 0);
	free(buffer);
	fprintf (stderr, "pucr: %u --> %u bytes\n",fileLen, outlen);
//	fwrite (out_buffer, sizeof(char), outlen, stdout);
	pucrunch_256_decode (output, &outlen, out_buffer, outlen);
	fwrite (output, sizeof(char), outlen, stdout);

	return(0);
}
