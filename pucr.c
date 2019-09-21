#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>


#include <time.h>
static clock_t start, end;
static double cpu_time_used;

#define LZ	0x0000
#define LIT	0x0002
#define RLE	0x0003
#define MARKED	0x0080

struct pucr_node {
	uint16_t	rle_count;
	uint16_t	lz_count;
	uint16_t	lz_offset;
	uint16_t	seen_before;
//	uint8_t		hash;

	uint16_t	way_to_go;
	uint16_t	next_node;
	uint32_t	bits_to_end;
	uint8_t 	new_esc;
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


static uint8_t int_log2 (uint32_t x) {
	if (x < 65535) {
		return (int_log2_lut[x]);
	} else {
		/* taken from Henry S. Warren Jr.'s 2nd edition of Hacker's Delight */
		uint8_t n = 0;
		if (x <= 0x0000FFFF) { n = n + 16; x = x << 16;}
		if (x <= 0x00FFFFFF) { n = n +  8; x = x <<  8;}
		if (x <= 0x0FFFFFFF) { n = n +  4; x = x <<  4;}
		if (x <= 0x3FFFFFFF) { n = n +  2; x = x <<  2;}
		n = n - (x >> 31);
		return (30 - n);
	}
}


// ---=== RLE RLE RLE RLE RLE RLE RLE RLE RLE RLE RLE RLE RLE RLE RLE RLE ===---

static uint32_t generate_rle_ranks (uint16_t rle_occ[], uint8_t rle_rank[], uint8_t *rank_len) {

	uint16_t max_index, max_value;
	int8_t i;

	for (i=0;i<15;i++) {
		/* find characters with highest number of occurences
		 * but not below 2 for ranks 3 and following */
		max_value = (i < 3)?1:2;
		max_index = 0xFFFF;
		for (int j=0; j < 256; j++) {
			if (rle_occ[j] > max_value) {
				/* check that candidate J is not already in list */
				int k;
				for (k=0;k<i;k++) if ( j==rle_rank[k] ) break;
				if (k == i) {
					max_value = rle_occ[j];
					max_index = j;
				}
			}
		}
		if ( max_index != 0xFFFF ) {
			rle_rank[i] = max_index;
			*rank_len = *rank_len + 1;
fprintf (stderr, "Rnk. %2u Occ. %5ux [%3u]\n",i, rle_occ[max_index], max_index ); // !!!
		} else {
			break;
		}
	}
	*rank_len = i;
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


static uint32_t find_matches (struct pucr_node graph[], const uint8_t *inbuf, uint16_t in_len) {

	int i,j,k;
	uint16_t lastocc[256][256];// NNN !!!

	for (i=0;i<256;i++) for (j=0;j<256;j++) lastocc[i][j] = 0xFFFF;

	for (i=0; i < in_len;i++ ) {
		graph[i].seen_before = lastocc[inbuf[i]][inbuf[i+1]];
		lastocc[inbuf[i]][inbuf[i+1]] = i;

// !!!	       graph[i].hash = inbuf[i]*3 + inbuf[i+1]*5 + inbuf[i]*7;

// we will be able to skip at least 2 characters when comparing
// because we know we already have a 2-byte match.
// even more if rle is longer in that position !!!
int off = (graph[i].rle_count -1);
off = off<=2?2:off;

		graph[i].lz_offset = 0xFFFF; // !!! neccessary for msb-lsb-offset calculation of number of bits
		graph[i].lz_count =  0x0000; // = no  match yet
		for (j=lastocc[inbuf[i]][inbuf[i+1]]; j != 0xFFFF; j=graph[j].seen_before) { // check at offset J for a match ...
			if (graph[i].rle_count != graph[j].rle_count) continue;
			k = j + off;
//if (graph[i+off-1].hash == graph[k-1].hash) {
			for (k=j+off; (k < i) && ((i+k-j)<in_len); k++) { // ... of length (K-J)
				if (inbuf[k] != inbuf[i+(k-j)]) break;
			}
//}
			if (k > i) k = j ; // no match
			k -= j; // match length
			if (k > 1) {
				// treatment of a match found
				// naive weighting function: the longer the better // replace with weighting/boundary conditions later on !!!
				if (k > graph[i].lz_count) {
					graph[i].lz_offset = i-j; // store the 'real' offset
					graph[i].lz_count=k;
				}
			}
		}
	}
}


static uint32_t find_optimal_lz_offset_no_lsb (struct pucr_node graph[], uint16_t in_len,
					       uint8_t no_esc,
					       uint8_t *no_lz_offset_lsb, uint8_t *no_lz2_offset_bits) {

	uint32_t lz_occ[16] = {0};
	uint32_t lz2_occ[16] = {0};
	// after all matches are finally determinded, calculate the optimum of lsb-msb coding ratio
	// first step: count them, especially the bits their offsets need
	for (int i=0; i != in_len; i=graph[i].next_node) {
		if (graph[i].way_to_go == LZ) {
			// lz_length is not necessarily the real way to go,
			// but next_node is
			if ((graph[i].next_node - i) >= 3) {
				// for the 3-byte and longer matches: calculate the case with "p" pure LSB and the rest E.Gamma-coded
			        for (int p=0; p <16; p++) // 0 = 1-bit offset, 1 = 2-bit offset, ...
					lz_occ[p] += (p+1) + int_log2( ((graph[i].lz_offset-(graph[i].next_node-i))>>(p+1))+1 )*2+1;
			} else {
				int p = int_log2 (graph[i].lz_offset - 2);
				p = p<=15?p:15;
				lz2_occ[p]++;
			}
		}
	}

	uint32_t min_sum = 0xFFFFFFFF;
	uint8_t min_no_lsb = 0;
	for (int m=15; m>=0; m--) {
		if (lz_occ[m] < min_sum) {
			min_no_lsb = m;
			min_sum = lz_occ[m];
		}
	}
	*no_lz_offset_lsb = min_no_lsb + 1;

	int max_p = 16 - 2 - no_esc - 1;
	int min_lz2_bits = int_log2 (in_len);
	min_lz2_bits = (min_lz2_bits <= 8)?min_lz2_bits:8;
	int max_saving = 0;
	for (int p=0; p<max_p; p++) {
		// count the bits saved for each LZ2 possible (i.e. if offset fits in p bits)
		int gained = 0;
		for (int q=0; q<=p; q++) gained += (lz2_occ[q] * (16 - (p+1 + no_esc +2)));
		// lost bits due to unused LZ2 opportunities as offset length is longer than p bits
		int lost = 0;
		for (int q=p+1; q<max_p; q++) lost += (lz2_occ[q] * (16 - (p+1 + no_esc +2)));
		// weigh against each other and store if minimum
		if ((gained-lost) > max_saving)  {
			max_saving = gained-lost;
			min_lz2_bits = p + 1;
		}
// fprintf (stderr, "LZ2 OFFSET %u BITS SAVED %d   LOST %d   RESULT: %d [maxp: %u]\n",p+1, gained, lost, gained-lost, max_p);
	}
// fprintf (stderr, "ESC: %2u LZ2 OFFSET %u BITS\n",no_esc, min_lz2_bits);
	*no_lz2_offset_bits = min_lz2_bits;
	return (0);
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
                              uint8_t no_esc, int *startEscape) {

    int i, j, states = (1<<no_esc);
    long minp = 0, minv = 0, other = 0;
    long a[256]; /* needs int/long */
    long b[256]; /* Remembers the # of escaped for each */
    int esc8 = 8-no_esc;

    for (i=0; i<256; i++)
        b[i] = a[i] = -1;

    if (states>256) {
        fprintf(stderr, "Escape optimize: only 256 states (%d)!\n",
                states);
        return 0;
    }

    for (i=0; i != in_len; i = graph[i].next_node)
	if (graph[i].way_to_go == LIT) graph[i].way_to_go = MARKED;

    for (i=in_len-1; i>=0; i--) {
        if (graph[i].way_to_go == MARKED) {
	    graph[i].way_to_go = LIT;
            int k = (inbuf[i] >> esc8);

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
    /* find smallest a[j] with smallest "smallest" escape code value */
    if (startEscape) {
        i = in_len;      /* make it big enough */
        for (j=states-1; j>=0; j--) {
            if (a[j] <= i) {
                *startEscape = (j << esc8);
                i = a[j];
            }
        }
    }
    return b[startEscape ? (*startEscape>>esc8) : 0] + 1;
}


static uint32_t find_best_path (struct pucr_node graph[], const uint8_t *inbuf, uint16_t in_len,
				uint16_t rle_occ[], uint8_t rle_rank[], uint8_t rank_len,
				uint8_t no_esc,
				uint8_t lz_lsb, uint8_t lz2_bits,
				uint8_t fast_lane) {

// !!! not for empty input -- sure ???
	/* create rle cost (per character) table w/o cost for following length value */
	uint32_t rle_cost_table[256];
	for (int j=0; j < 256; j++) {
		/* is it a ranked character? */
		int k;
		for (k=0; k < rank_len; k++) if ( j == rle_rank[k] ) break;
		if (k != rank_len) {
			// ranked
			rle_cost_table[j] = 3 + int_log2(k+1)*2+1; // + (32 / (3 * rle_occ[j] + 1)); // this is a try to take table costs into account !!!
		} else {
			// not ranked
			rle_cost_table[j] = 3 + int_log2(rank_len)+1 + 8;
		}
	}
	/* however, it is some blur left here: depending on choices in the following optimizer loop, rle occurences and their
	   cost may change again, which could affect optimizer choices again, which might change costs again, which might ... */

	// graph[in_len] is the (virtual) EOF symbol
	graph[in_len].bits_to_end = 0;
	/* check for each node the best way to go */
	for (int i=in_len-1; i >= 0; i--) {
		uint32_t lit_cost;
		uint32_t lz_cost, lz_way;
		uint32_t rle_cost, rle_way;
		uint32_t min = 0xFFFFFFFF;

		/* determine possible rle costs */
		/* this loops also checks for all shorter lengths */
		/* !!! ref org line 1044: only shorter ones that are a power of two? */
		for (int j=graph[i].rle_count; j > 1; j=fast_lane?1<<(int_log2(j)-1):j-1 ) {
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
		for (int j=graph[i].lz_count; j > 1; j=fast_lane?1<<(int_log2(j)-1):j-1) {
			// in case j == 2, the offset to be encoded may not be larger than 2^lz2_bits; otherwise lz is no option
			if ( (j == 2) && ((graph[i].lz_offset-2) >= (0x01 << lz2_bits)) ) break;
			lz_cost = (j==2)?2+lz2_bits:
				  int_log2(j-1)*2+1 + lz_lsb + int_log2(((graph[i].lz_offset-j)>>lz_lsb)+1)*2+1;
			lz_cost += graph[i+j].bits_to_end;
			if (lz_cost < min) {
				min = lz_cost;
				lz_way = i+j;
			}
		}
		// only if a minimum was found in the loop, add the number of ESCape bits
		// otherwise, keep up the high cost of 0xFFFFFFFF
		// !!! (min+1) == 0 is a bit hacky and relies on data type properties
		lz_cost = (min+1)?min+no_esc:min;

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


static int outPointer = 0;
static uint8_t bitMask = 0x80;
static int bitCount = 0;
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


static void put_value (int value) {
    int bits = 0, count = 0;

    while (value>1) {
        bits = (bits<<1) | (value & 1); /* is reversed compared to value */
        value >>= 1;
        count++;
        put_bit (1);
    }
// !!!    if (count<maxGamma) // no max_gamma yet !!!
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


static uint32_t output_path (struct pucr_node graph[], const uint8_t *inbuf, uint16_t in_len,
			uint16_t rle_occ[], uint8_t rle_rank[], uint8_t rank_len,
			uint8_t start_esc, uint8_t no_esc,
			uint8_t lz_lsb, uint8_t lz2_bits) {

	uint8_t graph_current_esc = start_esc;
	uint8_t current_esc = start_esc;
	uint8_t current_esc_mask = 0xFF00 >> no_esc;
	uint8_t current_esc_8 = current_esc >> (8 - no_esc);

	// output header
	// decompressed length
	put_value ((in_len>>8)+1);
	put_n_bits (in_len,8);
	// number of ESCape bits and first ESCape sequence
	put_n_bits (no_esc,4);
	put_n_bits (current_esc_8,no_esc);
	// RLE table count and table w/ up to 15 entries
	put_n_bits (rank_len,4);
	for (uint8_t i=0; i < rank_len; i++)
		put_n_bits(rle_rank[i],8);
	// number of bits for LZ2 offset
	put_n_bits (lz2_bits,4);
	// number of LSBs for LZ offset
	put_n_bits (lz_lsb,4);

	// output data
	uint8_t no_rle_char_esc_bits = (rank_len==0)?0 : int_log2(rank_len)+1;
fprintf (stderr, "RLE CHAR ESC CODE: no of 1's: %u\n", no_rle_char_esc_bits);

	for (int i=0; i!=in_len; i=graph[i].next_node) {

		switch (graph[i].way_to_go) {

			case RLE:
//				fprintf (stderr, "%u. RLE(%c,%u) --> %u [%u]",i,inbuf[i], graph[i].next_node - i, graph[i].next_node,  graph[i].bits_to_end - graph[graph[i].next_node].bits_to_end);
				// ESC ++ 011
				put_n_bits (current_esc_8, no_esc);
				put_bit (0); put_bit (1); put_bit (1);
                                // output RLE count (real count minus 1, i.e. 1 instead of 2 etc.) E. Gamma coded
                                put_value (graph[i].next_node - i - 1);
				// check if it is a ranked character
				int k;
				for (k=0; k < rank_len; k++) if ( inbuf[i] == rle_rank[k] ) break;
				if (k != rank_len) {
					// is ranked: E. Gamma coded pointer to rank list
					put_value (k+1);
//					fprintf (stderr, " {%u} ", no_esc + 3 + int_log2(k+1)*2+1 + int_log2(graph[i].next_node - i - 1)*2+1);
				} else {
					// not ranked: determine '11..11'-code length (already determined outside loop) followed by character
					put_n_bits (0xFF, no_rle_char_esc_bits);
					put_n_bits (inbuf[i], 8);
//					fprintf (stderr, " {%u} ", no_esc + 3 + no_rle_char_esc_bits + 8 + int_log2(graph[i].next_node - i - 1)*2+1);
				}
				break;
			case LZ:
/*				fprintf (stderr, "%u. LZ(%u,%u) --> %u [%u]",i,
										graph[i].next_node-i,
										graph[i].lz_offset,
										graph[i].next_node,
										-graph[graph[i].next_node].bits_to_end + graph[i].bits_to_end);
*/				put_n_bits (current_esc_8, no_esc);
				if ( (graph[i].next_node-i) == 2 ) {
					// LZ2: ESC + 000 + LZ2 OFFSET (lz2_bits)
	                                put_bit (0); put_bit (0); // put_bit (0);
					put_n_bits (graph[i].lz_offset-2, lz2_bits);
//					fprintf (stderr, " {%u} ", no_esc + 2 + lz2_bits);
				} else {
					// LZ: ESC + E. Gamma coded actually used lz length (next_node - i - 1), not neccessariliy lz_len ...
					put_value (graph[i].next_node - i - 1);
					// ... and also E. Gamma coded LZ offset minus the actually used length of lz match(=next_node-i)
					int off = (graph[i].lz_offset-(graph[i].next_node-i));
					// MSBs (+1 for E. Gamma code)
					put_value ((off >> lz_lsb)+1);
					// LSBs
					put_n_bits (off, lz_lsb);
//					fprintf (stderr, " {%u} ", no_esc + int_log2(graph[i].next_node-i-1)*2+1 + int_log2((off>>lz_lsb)+1)*2+1 + lz_lsb);
				}
				break;
			case LIT:
//				fprintf (stderr, "%u. LIT(%c) --> %u [%u]",i,inbuf[i],graph[i].next_node, -graph[graph[i].next_node].bits_to_end + graph[i].bits_to_end+(graph_current_esc == graph[i].new_esc?0:no_esc+3));
				if (graph[i].new_esc != graph_current_esc) {
//					fprintf (stderr, "  NEW possible ESC: %u", graph[i].new_esc);
					graph_current_esc = graph[i].new_esc;
				}
				//
				if ( (inbuf[i] & current_esc_mask) != current_esc) {
					put_n_bits (inbuf[i], 8);
//					fprintf (stderr, " {%u} ",8);
				} else {
					// ESC = first bits of LITeral
					put_n_bits (current_esc_8, no_esc);
					// 010
	                                put_bit (0); put_bit (1); put_bit (0);
					// new ESC code
					current_esc = graph[i].new_esc;
					current_esc_8 = current_esc >> (8 - no_esc);
					put_n_bits (current_esc_8, no_esc);
//					fprintf (stderr, "NEW ESCAPE: %u ", current_esc);
					// rest of LITeral
					put_n_bits (inbuf[i], 8- no_esc);
//					fprintf (stderr, " {%u} ",8+3+no_esc);
				}
				break;
			default:
				// !!! shall not happen
				break;
		}
//		printf ("\n");
	}
//	printf ("OUT POINTER: %u  BITMASK: %u\n", outPointer, bitMask);
	flush_bits();
//	printf ("OUT POINTER: %u\n", outPointer);
}


uint32_t pucrunch_256_encode (uint8_t *outbuf,uint32_t *out_len,
                             const uint8_t *inbuf, uint16_t in_len,
			     uint8_t fast_lane) {

start_clock();

	struct pucr_node graph[65536];
	uint16_t rle_occ[256];
	uint8_t rle_rank[15];
	uint8_t rank_len;

	// reset_graph_path (graph, in_len);
	count_rle_and_generate_occ_from_inbuf (graph, inbuf, in_len, rle_occ);
	generate_rle_ranks (rle_occ, rle_rank, &rank_len);

fprintf (stderr, "--- COUNTING RLE: ");
print_clock(); start_clock();
	find_matches (graph, inbuf, in_len); // !!! option for different weighting
fprintf (stderr, "--- FINDING MATCHES: ");
print_clock(); start_clock();

	uint8_t lz_lsb = int_log2 (in_len)*2/3;
	uint8_t lz2_bits = int_log2 (in_len)*2/3;
		lz2_bits = lz2_bits>8?8:lz2_bits;
//	uint8_t lz2_bits = 8;

	uint32_t min_no_bits = 0xFFFFFFFF;
	uint8_t min_no_esc;
	uint8_t min_lz_lsb;
	uint8_t min_lz2_bits;
	int min_start_esc;

	uint8_t no_esc;
	int start_esc = 0;

	if (fast_lane < 3) {
		for (no_esc = (fast_lane<2)?0:1; no_esc < ((fast_lane<2)?9:4); no_esc++) {
			// the first path generated using sane default
			find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rank_len, no_esc, lz_lsb, lz2_bits, fast_lane);
			find_optimal_lz_offset_no_lsb (graph, in_len, no_esc, &lz_lsb, &lz2_bits);
//			find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rank_len, no_esc, lz_lsb, lz2_bits, fast_lane);
			int escaped = optimize_esc (graph, inbuf, in_len, no_esc, &start_esc);

if (fast_lane < 2)		update_rle_occ_from_graph (graph, inbuf, in_len, rle_occ);
if (fast_lane < 2)		generate_rle_ranks (rle_occ, rle_rank, &rank_len);
if (fast_lane < 2)		find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rank_len, no_esc, lz_lsb, lz2_bits, fast_lane);

	            	/* Compare value: bits lost for escaping -- bits lost for prefix */
			if ((graph[0].bits_to_end + (no_esc+3) * escaped) < min_no_bits) {
				min_no_esc = no_esc;
				min_no_bits = graph[0].bits_to_end + (no_esc+3)*escaped; // "- escaped" is for those counted 9-bit LIT in find_best_path
				// also store the coresponding lz_lsb value
				min_lz_lsb = lz_lsb;
				min_lz2_bits = lz2_bits;
				min_start_esc = start_esc & ((0xff00>>no_esc) & 0xff);
			}
fprintf (stderr, "--- OPTIMIZING ESC %u: ", no_esc);
print_clock(); start_clock();
		}
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
	find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rank_len, min_no_esc, min_lz_lsb, min_lz2_bits, fast_lane);
if (fast_lane < 2)		update_rle_occ_from_graph (graph, inbuf, in_len, rle_occ);
if (fast_lane < 2)		generate_rle_ranks (rle_occ, rle_rank, &rank_len);
if (fast_lane < 2)		find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rank_len, min_no_esc, min_lz_lsb, min_lz2_bits, fast_lane);
	int escaped = optimize_esc (graph, inbuf, in_len, min_no_esc, &min_start_esc);
	min_no_bits = graph[0].bits_to_end + (min_no_esc+3)*escaped; // "- escaped" is for those counted 9-bit LIT in find_best_path

	output_path (graph, inbuf, in_len, rle_occ, rle_rank, rank_len, min_start_esc, min_no_esc, min_lz_lsb, min_lz2_bits );

fprintf(stderr, "-----== encoding stats ==-----\n#ESC: 			%u\n#ESCAPED LIT:		%u\n#LZ OFFSET LSB:		%u\n#LZ2 OFFSET BITS	%u\n#RANKED RLEs		%u\nIN:				%u\nOUT w/o header:		%u\n------------------------------\n",
	min_no_esc, escaped, min_lz_lsb, min_lz2_bits, rank_len, in_len, (min_no_bits + 7)>>3);

	*out_len = (min_no_bits + 7) >> 3;
fprintf (stderr, "--- GENERATING FINAL PATH: ");
print_clock(); start_clock();

	return (0);
}

// ---=== DECODE DECODE DECODE DECODE DECODE ===---

static const unsigned char *up_Data;
static int up_Mask, up_Byte;


static void up_SetInput(const unsigned char *data) {
    up_Data = data;
    up_Mask = 0x80;
    up_Byte = 0;
}


static int up_GetBits(int bits) {
    int val = 0;

    while (bits--) {
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

    while (1) { // !!! was: (i<maxGamma) {
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


static uint32_t pucrunch_256_decode (uint8_t *outbuf,uint32_t *out_len,
                             const uint8_t *inbuf, uint32_t in_len) {

start_clock();
	up_SetInput (inbuf);

	// decompressed length
	// datatype limits block size, needs to be signed for loop
	int32_t len = (up_GetValue()-1)<<8;
	len += up_GetBits(8);
	// number of ESCape bits and first ESCape sequence
	uint8_t no_esc = up_GetBits(4);
	uint8_t current_esc8 = up_GetBits(no_esc);
	// RLE table count and table w/ up to 15 entries
	uint8_t rank_len = up_GetBits(4);
	uint8_t rle_rank[15];
	for (uint8_t i=0; i < rank_len; i++)
		rle_rank[i] = up_GetBits(8);
	// number of '1' bits to identify a following character (as opposed to rank pointer)
	uint8_t no_rle_char_esc_bits = (rank_len==0)?0 : int_log2(rank_len)+1;

	// number of bits for LZ2 offset
	uint8_t lz2_bits = up_GetBits(4);
	// number of LSBs for LZ offset
	uint8_t lz_lsb = up_GetBits(4);

fprintf(stderr, "-----== decoding stats ==-----\n#ESC: 			%6u\n#LZ OFFSET LSB:		%6u\n#LZ2 OFFSET BITS	%6u\n#RANKED RLEs		%6u\nIN:			%6u\nOUT w/o header:		%6u\n------------------------------\n",
	no_esc, lz_lsb, lz2_bits, rank_len, in_len, len);

	uint8_t current_esc = current_esc8 << (8 - no_esc);
	uint8_t esc_mask = 0xFF00 >> no_esc;
uint16_t rle_ranked_cnt = 0;
	for (uint16_t cnt=0;cnt<len;) {
		uint8_t first_fetch = up_GetBits(no_esc);
		if (first_fetch == current_esc8) {
			uint16_t fetch = up_GetValue();
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
						// ESCaped LITeral
						// get the new ESCape code
						current_esc8 = up_GetBits(no_esc);
						current_esc = current_esc8 << (8 - no_esc);
						// fetch the rest of the LITeral
						first_fetch = (first_fetch<<(8-no_esc)) | up_GetBits(8-no_esc);
						outbuf[cnt++] = first_fetch;
//fprintf (stderr, "%u. LIT (%c)  --ESCAPED--\n", cnt-1, first_fetch);
					}
				} else {
					// LZ - to be specific, this is LZ2
					uint16_t offset = up_GetBits (lz2_bits);
// fprintf (stderr, "%u. LZ (2,%u)\n", cnt, offset+2);
					outbuf[cnt]=outbuf[cnt-2-offset]; cnt++;
					outbuf[cnt]=outbuf[cnt-2-offset]; cnt++;
				}
			} else {
				// common LZ, the fetch already contains the LZ count
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
			first_fetch = (first_fetch<<(8-no_esc)) | up_GetBits(8-no_esc);
			outbuf[cnt++] = first_fetch;
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
	uint32_t outlen = 65535;

	pucrunch_256_encode (output, &outlen, buffer, fileLen, 0);
	free(buffer);
	fprintf (stderr, "pucr: %u --> %u bytes\n",fileLen, outlen);

//	fwrite (out_buffer, sizeof(char), outlen, stdout);
	pucrunch_256_decode (output, &outlen, out_buffer, outlen);
	fwrite (output, sizeof(char), outlen, stdout);

	return(0);
}
