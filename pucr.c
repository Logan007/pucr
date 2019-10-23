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
	uint16_t	lz_max_offset;
	uint16_t	lz_cur_offset;
	uint16_t	seen_before;
//	uint8_t		hash;

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


// ---=== MOVE TO FRONT   MOVE TO FRONT   MOVE TO FRONT   MOVE TO FRONT ===---

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


// ---=== HUFFMAN HUFFMAN HUFFMAN HUFFMAN HUFFMAN HUFFMAN HUFFMAN HUFFMAN  ===---

struct huffman_code_t {
    uint8_t code[32];
    uint8_t code_len;
};

struct huffman_node {
    int32_t count, l, r;
    bool used;
    struct huffman_code_t code;
};


uint32_t huffman_tree_encode (struct huffman_node table[], uint16_t current_node,
                              struct huffman_code_t current_code,
                              uint8_t *encoded_tree, uint16_t *tree_position,
                              uint8_t no_symbol_bits) {

    if (table[current_node].l == -1) {
        // ouput '1' as code for a following leaf
        encoded_tree[*tree_position >> 3] = encoded_tree[*tree_position >> 3] | (0x80 >> (*tree_position & 0x07));
        *tree_position = *tree_position + 1;
        uint8_t cn = current_node;
        // write 'no_symbol_bits'-bit uncompressed value
        for (int8_t i=no_symbol_bits-1; i >= 0; i--) {
            if ((cn >> i) & 0x01) {
	        encoded_tree[*tree_position >> 3] = encoded_tree[*tree_position >> 3] | (0x80 >> (*tree_position & 0x07));
	    } else {
	        encoded_tree[*tree_position >> 3] = encoded_tree[*tree_position >> 3] & (0xFF00 >> (*tree_position & 0x07));
	    }
            *tree_position = *tree_position + 1;
        }

        // complete missing shifts and write code to leaf
        for (uint16_t i=current_code.code_len; (i & 0x07) != 0; i++)
            current_code.code[i >> 3] = current_code.code[i >> 3] << 1;
        for (uint8_t i=0; i < 32; i++)
            table[current_node].code.code[i] = current_code.code[i];
        table[current_node].code.code_len = current_code.code_len;
    } else {
        // output'0' to indicate following l and r which are the recursively determined respective sub-trees (0lr)
        encoded_tree[*tree_position >> 3] = encoded_tree[*tree_position >> 3] & (0xFF00 >> (*tree_position & 0x07));
        *tree_position = *tree_position + 1;

        // construct codeword (l=0)
        current_code.code[current_code.code_len >> 3] = current_code.code[current_code.code_len >> 3] << 1;
        current_code.code_len++;
        // recursively parse left sub-tree
        huffman_tree_encode (table, table[current_node].l, current_code, encoded_tree, tree_position, no_symbol_bits);

        // construct codeword (r=1)
        current_code.code_len--;
        current_code.code[current_code.code_len >> 3] = current_code.code[current_code.code_len >> 3] | 0x01;
        current_code.code_len++;
        // recursively parse right sub-tree
        huffman_tree_encode (table, table[current_node].r, current_code, encoded_tree, tree_position, no_symbol_bits);
    }
}


uint32_t generate_one_huffman_tree_and_part_of_codebook (uint8_t huffman_tree[], uint16_t *huffman_tree_len_bits,
                                                         struct huffman_node table[], uint16_t *first_head,
                                                         uint8_t no_symbol_bits, uint8_t preceding_bits) {


// !!! max size of encoded huffman tree: 320 byte ; do we need to check for it?
// !!! max size of huffman_node table: 256 leafs + max 255 nodes = 511 entries
// !!! head on first call = 255

        /* the huffman table comes prefilled with bytes' occurences in table[i].count */
        /* mark those entries without hits as 'already used' so they do not get checked later */
        // also: those that do not match the leading mask
        uint8_t preceding_bits8 = preceding_bits << no_symbol_bits;
        uint8_t mask = 0xFF00 >> (8 - no_symbol_bits);

        for (uint16_t i=0; i < 256; i++) {
            if ( (table[i].count == 0) | ((i & mask) != preceding_bits8) )
                table[i].used = true;
            else
                table[i].used = false;
            table[i].l = -1;
            table[i].r = -1;
        }
        for (uint16_t i=256; i<511; i++) table[i].used=false;

        uint16_t head = *first_head;
        do {
            head++;

            // LEFT SIDE
            uint32_t min = 0xFFFFFFFF;
            uint32_t min_idx = -1;
	    // skip nodes from formerly created trees
            for (uint16_t i=0; i < head; i=(i==255)?*first_head+1:i+1) {
                if (!table[i].used) {
                    if (table[i].count < min) {
                        min = table[i].count;
                        min_idx = i;
                    }
                }
            }

            table[head].l = min_idx;

	    // not even a single node (empty tree)
	    if (!(min+1))
		break;

            table[head].count = table[min_idx].count;
            table[min_idx].used = true;

	    // RIGHT SIDE
            min = 0xFFFFFFFF;
            min_idx = -1;
            // skip nodes from formerly created trees
            for (uint16_t i=0; i < head; i=(i==255)?*first_head+1:i+1) {
                if (!table[i].used) {
                    if (table[i].count < min) {
                        min = table[i].count;
                        min_idx = i;
                    }
                }
            }
            if (min+1) {
                table[head].r = min_idx;
                table[head].count += table[min_idx].count;
                table[min_idx].used = true;
            } else {
                head=table[head].l;
                break;
            }
        } while (true);

        *first_head = head<256?*first_head:head;

        /* output tree and generate codebook */
        struct huffman_code_t empty_code = {};
	empty_code.code_len = 0;
//	empty_code.code [0] = 0;
        huffman_tree_encode (table, head, empty_code, huffman_tree, huffman_tree_len_bits, no_symbol_bits);
}


uint32_t generate_huffman_trees_and_codebook (struct pucr_node graph[], uint8_t const *inbuf, uint16_t in_len, // !!! was const !!! inbuf
                                              uint8_t huffman_tree[], uint16_t *huffman_tree_len_bits,
                                              struct huffman_node table[], uint16_t head_array[],
                                              uint8_t no_symbol_bits) {

	uint32_t ret = 0;
        uint8_t esc_mask = 0xFF00 >> (8-no_symbol_bits);
        uint8_t mask = 0x00FF >> (8-no_symbol_bits);

// ---

//        for (uint16_t i=0; i < 256; i++) table[i].count = 0;
//        // count LITeral occurences
//        for (uint16_t i=0; i < in_len; i=graph[i].next_node)
//                if (graph[i].way_to_go == LIT)
//                        table[  graph[i].cur_lit&mask  ].count++;
//        *huffman_tree_len_bits = 0; uint16_t head = 255;
//
//        generate_one_huffman_tree_and_part_of_codebook (huffman_tree, huffman_tree_len_bits,
//                                                         table, &head,
//                                                         8, 0);
//
//        int32_t sum = 0;
//        for (uint16_t j=0; j < mask+1; j++)
//                 sum += table[j].count * (no_symbol_bits - table[j].code.code_len);
//        fprintf (stderr, "A SINGLE (!) HUFFMAN FOR ALL SYMBOLS WITH %u BIT SYMBOLS SAVED %d BITS\n",no_symbol_bits, sum - *huffman_tree_len_bits);
//        fprintf (stderr, "HEAD: %u    TREE LEN: %u BITS\n",head,*huffman_tree_len_bits);
// ------

        for (uint16_t i=0; i < 256; i++) table[i].count = 0;
        // count LITeral occurences
        for (uint16_t i=0; i < in_len; i=graph[i].next_node)
                 if (graph[i].way_to_go == LIT)
                        table[graph[i].cur_lit].count++;

	uint16_t head = 255;

	for (uint16_t i=0; i<320;i++) huffman_tree[i]=0;

//	uint8_t huffman_tree[320];
//       uint16_t huffman_tree_len_bits = 0;
//        *huffman_tree_len_bits = 0;


// is 320 sufficent (including the indicators?) !!! ??? it should

        for (uint16_t i=0; i<(uint8_t)(1<<(8-no_symbol_bits)); i++) {
                uint16_t old_tree_len_bits = *huffman_tree_len_bits;
                // leave placeholder space for a 1-bit indicator
		*huffman_tree_len_bits = *huffman_tree_len_bits + 1;
		generate_one_huffman_tree_and_part_of_codebook (huffman_tree, huffman_tree_len_bits,
                                                         table, &head,
                                                         no_symbol_bits, i);
                int32_t sum = 0;
                for (uint16_t j=(i << no_symbol_bits); j < ((i+1) << no_symbol_bits); j++)
                         sum += table[j].count * (no_symbol_bits - table[j].code.code_len);
// fprintf (stderr, "--- HUFFMAN WITH %u BIT SYMBOLS WOULD SAVE %d BITS\n",no_symbol_bits, sum-(*huffman_tree_len_bits-old_tree_len_bits));
		// if positive savings... (even the 1-bit indicator needs to pay off)
		if ( (*huffman_tree_len_bits-old_tree_len_bits) < (sum+1) ) {
			// set indicator to '1'
			huffman_tree [old_tree_len_bits >> 3] = huffman_tree [old_tree_len_bits >> 3] | (0x80 >> (old_tree_len_bits & 0x07));
			// store the head for easier access to the table while encoding
			head_array[i] = head;
			ret += sum-(*huffman_tree_len_bits-old_tree_len_bits);
		} else {
			// set indicator to '0' (do not care about the following bits)
			huffman_tree [old_tree_len_bits >> 3] = huffman_tree [old_tree_len_bits >> 3] & (0xFF7F >> (old_tree_len_bits & 0x07));
			// do not use the tree for this part (but leave the '0' bit indicator)
 			*huffman_tree_len_bits = old_tree_len_bits +1;
			// indicate 'no tree' to easen LITerals' encoding
			head_array[i] = 0xFFFF;
		}
        }
	return (ret);
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
// fprintf (stderr, "Rnk. %2u Occ. %5ux [%3u]\n",i, rle_occ[max_index], max_index ); // !!!
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


// ---=== LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ LZ ===---


static uint32_t find_closer_match (struct pucr_node graph[], const uint8_t *inbuf, uint16_t in_len,
				   uint16_t pattern_start, uint16_t pattern_length, uint16_t min_match_start,
                            	   uint16_t last_occ[256][256]) {

        int j,k ;

	// we will be able to skip at least 2 characters when comparing
	// because we know we already have a 2-byte match.
	// even more if rle is longer in that position !!!
	int off = (graph[pattern_start].rle_count -1);
	off = off<=2?2:off;
//	graph[pattern_start].lz_cur_offset = graph[pattern_start].lz_max_offset; // 0xFFFF;
        for (j=last_occ[inbuf[pattern_start]][inbuf[pattern_start+1]]; (j >= min_match_start) && (k != 0xFFFF); j=graph[j].seen_before) { // check at offset J for a match ...
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


static uint32_t find_matches (struct pucr_node graph[], const uint8_t *inbuf, uint16_t in_len, uint16_t last_occ[256][256]) {

	int i,j,k;

	for (i=0;i<256;i++) for (j=0;j<256;j++) last_occ[i][j] = 0xFFFF;

	for (i=0; i < in_len;i++ ) {
		graph[i].seen_before = last_occ[inbuf[i]][inbuf[i+1]];
		last_occ[inbuf[i]][inbuf[i+1]] = i;

// !!!	       graph[i].hash = inbuf[i]*3 + inbuf[i+1]*5 + inbuf[i]*7;

// we will be able to skip at least 2 characters when comparing
// because we know we already have a 2-byte match.
// even more if rle is longer in that position !!!
int off = (graph[i].rle_count -1);
off = off<=2?2:off;

		graph[i].lz_max_offset = 0xFFFF; // !!! neccessary for msb-lsb-offset calculation of number of bits
		graph[i].lz_count =  0x0000; // = no  match yet
		for (j=last_occ[inbuf[i]][inbuf[i+1]]; j != 0xFFFF; j=graph[j].seen_before) { // check at offset J for a match ...
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
					graph[i].lz_max_offset = i-j; // store the 'real' offset
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
					lz_occ[p] += (p+1) + int_log2( ((graph[i].lz_cur_offset-(graph[i].next_node-i))>>(p+1))+1 )*2+1;
			} else {
				int p = int_log2 (graph[i].lz_cur_offset - 2);
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
                              uint8_t no_esc, int *startEscape,
			      uint8_t mtf, uint8_t mtf_second_line) {

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

    uint8_t alphabet[256]; for (i=0; i<256;i++) alphabet[i]=i;

    for (i=0; i != in_len; i = graph[i].next_node)
	if (graph[i].way_to_go == LIT) {
	    graph[i].way_to_go = MARKED;
	    // store the (maybe to be mtf-encoded) literal in the graph to preserve inbuf
	    graph[i].cur_lit = mtf?move_to_front_encode (inbuf[i], alphabet, mtf_second_line):inbuf[i];
	}

    for (i=in_len-1; i>=0; i--) {
        if (graph[i].way_to_go == MARKED) {
	    graph[i].way_to_go = LIT;
	    // use mtf-encoded literal for ESCape code evaluation
            int k = (graph[i].cur_lit >> esc8);

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
				uint8_t lz_lsb, uint8_t lz2_bits, uint16_t last_occ[256][256],
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
			if (rank_len == 0) rle_cost_table[j]--;
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
			if (fast_lane == 0) {
				int used_match_length = lz_way-i;
				if (used_match_length<graph[i].lz_count) { // LZ2 matches only in fast_lane 0
					find_closer_match (graph, inbuf, in_len, i, used_match_length,i-graph[i].lz_max_offset+1, last_occ);
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
			uint16_t *out_len,
			uint16_t rle_occ[], uint8_t rle_rank[], uint8_t rank_len,
			uint8_t start_esc, uint8_t no_esc, uint8_t mtf, uint8_t mtf_second_line,
			uint8_t lz_lsb, uint8_t lz2_bits,
			struct huffman_node table[], uint16_t head_array[], uint8_t huffman_tree[], uint16_t huffman_tree_len_bits) {

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
	// move to front
	if (mtf) put_bit (1); else put_bit (0);
	if (mtf) put_n_bits (mtf_second_line,8);
	// huffman tree(s)
	uint16_t i;
	for (i=0; i < huffman_tree_len_bits; i=i+8) {
		uint8_t no_bits = (uint16_t)(huffman_tree_len_bits-i)>=8?8:(huffman_tree_len_bits-i);
		put_n_bits (huffman_tree[i >> 3]>>(8-no_bits), no_bits );
	}

	// derive further parameters
	uint8_t no_rle_char_esc_bits = (rank_len==0)?0 : int_log2(rank_len)+1;
// fprintf (stderr, "RLE CHAR ESC CODE: no of 1's: %u\n", no_rle_char_esc_bits);

	// output data
	for (int i=0; i!=in_len; i=graph[i].next_node) {
// int old = bitCount;
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
//				fprintf (stderr, "%u. LZ(%u,%u) --> %u [%u]",i,
//										graph[i].next_node-i,
//										graph[i].lz_max_offset,
//										graph[i].next_node,
//										-graph[graph[i].next_node].bits_to_end + graph[i].bits_to_end);
				put_n_bits (current_esc_8, no_esc);
				if ( (graph[i].next_node-i) == 2 ) {
					// LZ2: ESC + 000 + LZ2 OFFSET (lz2_bits)
	                                put_bit (0); put_bit (0); // put_bit (0);
					put_n_bits (graph[i].lz_cur_offset-2, lz2_bits);
//					fprintf (stderr, " {%u} ", no_esc + 2 + lz2_bits);
				} else {
					// LZ: ESC + E. Gamma coded actually used lz length (next_node - i - 1), not neccessariliy lz_len ...
					put_value (graph[i].next_node - i - 1);
					// ... and also E. Gamma coded LZ offset minus the actually used length of lz match(=next_node-i)
					int off = (graph[i].lz_cur_offset-(graph[i].next_node-i));
					// MSBs (+1 for E. Gamma code)
					put_value ((off >> lz_lsb)+1);
					// LSBs
					put_n_bits (off, lz_lsb);
//					fprintf (stderr, " {%u} ", no_esc + int_log2(graph[i].next_node-i-1)*2+1 + int_log2((off>>lz_lsb)+1)*2+1 + lz_lsb);
				}
/* fprintf (stderr, "%u. LZ(%u,%u) --> %u [%u vs bitCount %u]\n",i,
										graph[i].next_node-i,
										graph[i].lz_max_offset,
										graph[i].next_node,
										-graph[graph[i].next_node].bits_to_end + graph[i].bits_to_end, bitCount-old);
*/				break;
			case LIT:
//				fprintf (stderr, "%u. LIT(%c) --> %u [%u]",i,inbuf[i],graph[i].next_node, -graph[graph[i].next_node].bits_to_end + graph[i].bits_to_end+(graph_current_esc == graph[i].new_esc?0:no_esc+3));
				if (graph[i].new_esc != graph_current_esc) {
//					fprintf (stderr, "  NEW possible ESC: %u", graph[i].new_esc);
					graph_current_esc = graph[i].new_esc;
				}
				if ( (graph[i].cur_lit & current_esc_mask) != current_esc) {
					// output the mtf-encoded literal
// w/o HUFFMAN:				put_n_bits (graph[i].cur_lit, 8);
					put_n_bits (graph[i].cur_lit >> (8-no_esc), no_esc);
					// select huffman tree if available
					// eventually, use codebook for output
					if ( head_array[(graph[i].cur_lit) >>(8-no_esc)] == 0xFFFF)
						put_n_bits (graph[i].cur_lit, 8- no_esc);
					else {
						for (int j=0; j < table[graph[i].cur_lit].code.code_len; j++) {
							if ( table[graph[i].cur_lit].code.code[j >> 3] & (0x80 >> (j & 7)) )
								put_bit(1); 
							else 
								put_bit(0);

						}

					}
// fprintf (stderr, " {%u} ",8);
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
					// rest of mtf-encoded LITeral:
					// select huffman tree if available
					// eventually, use codebook for output

					if ( head_array[(graph[i].cur_lit) >>(8-no_esc)] == 0xFFFF)
						put_n_bits (graph[i].cur_lit, 8- no_esc);
					else {
	                                        for (int j=0; j < table[graph[i].cur_lit].code.code_len; j++) {
	                                                if ( table[graph[i].cur_lit].code.code[j >> 3] & (0x80 >> (j & 7)) )
	                                                        put_bit(1); else put_bit(0);
	                                        }
					}
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
}


uint32_t pucrunch_256_encode (uint8_t *outbuf,uint16_t *out_len,
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

	uint16_t last_occ[256][256];
	find_matches (graph, inbuf, in_len, last_occ); // !!! option for different weighting

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

	// do the graph building with (standard) mtf turned on
	uint8_t mtf = 1;

	uint8_t no_esc;
	int start_esc = 0;

	if (fast_lane < 3) {
		for (no_esc = (fast_lane<2)?0:1; no_esc < ((fast_lane<2)?9:4); no_esc++) {
			// the first path generated using sane default
			find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rank_len, no_esc, lz_lsb, lz2_bits, last_occ, fast_lane);
			find_optimal_lz_offset_no_lsb (graph, in_len, no_esc, &lz_lsb, &lz2_bits);
//			find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rank_len, no_esc, lz_lsb, lz2_bits, last_occ, fast_lane);
			int escaped = optimize_esc (graph, inbuf, in_len, no_esc, &start_esc, mtf, 0);

if (fast_lane < 2)	update_rle_occ_from_graph (graph, inbuf, in_len, rle_occ);
if (fast_lane < 2)	generate_rle_ranks (rle_occ, rle_rank, &rank_len);
if (fast_lane < 2)	find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rank_len, no_esc, lz_lsb, lz2_bits, last_occ, fast_lane);

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
	find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rank_len, min_no_esc, min_lz_lsb, min_lz2_bits, last_occ, fast_lane);
if (fast_lane < 2)		update_rle_occ_from_graph (graph, inbuf, in_len, rle_occ);
if (fast_lane < 2)		generate_rle_ranks (rle_occ, rle_rank, &rank_len);
if (fast_lane < 2)		find_best_path (graph, inbuf, in_len, rle_occ, rle_rank, rank_len, min_no_esc, min_lz_lsb, min_lz2_bits, last_occ, fast_lane);

	// handling of remaining LITerals:
	int escaped;
	int huffed;
        uint8_t huffman_tree[320];
        uint16_t huffman_tree_len_bits = 0;
        struct huffman_node table[511];
	uint16_t head_array[256];

fprintf (stderr, "GENERATING FINAL GRAPH: ");
print_clock(); start_clock();

	// determine if we would be better off with or without mtf
	uint8_t min_mtf = 0;
	uint32_t min_mtf_bits = 0xFFFFFFFF;
	for (mtf = 0; mtf < 2; mtf++) {
		escaped = optimize_esc (graph, inbuf, in_len, min_no_esc, &min_start_esc, mtf, 0);
		min_no_bits = graph[0].bits_to_end + (min_no_esc+3)*escaped; // "- escaped" is for those counted 9-bit LIT in find_best_path
	        huffman_tree_len_bits = 0;
	        huffed = generate_huffman_trees_and_codebook (graph, inbuf, in_len, huffman_tree, &huffman_tree_len_bits, table, head_array, 8-min_no_esc);
		if ( (min_no_bits-huffed) < min_mtf_bits ) {
			min_mtf = mtf;
			min_mtf_bits = (min_no_bits - huffed);
		}
	}

	// if mtf, will a higher mtf_second_line be beneficial?
	uint8_t min_mtf_second_line = 0;
	uint32_t min_mtf_second_line_bits = 0xFFFFFFFF;
	if (min_mtf) {
		for (uint8_t mtf_second_line=0; mtf_second_line < 254; mtf_second_line++) {
			escaped = optimize_esc (graph, inbuf, in_len, min_no_esc, &min_start_esc,min_mtf , mtf_second_line);
			min_no_bits = graph[0].bits_to_end + (min_no_esc+3)*escaped;
		        huffman_tree_len_bits = 0;
		        huffed = generate_huffman_trees_and_codebook (graph, inbuf, in_len, huffman_tree, &huffman_tree_len_bits, table, head_array, 8-min_no_esc);
			if ( (min_no_bits-huffed) < min_mtf_second_line_bits ) {
				min_mtf_second_line = mtf_second_line;
				min_mtf_second_line_bits = (min_no_bits - huffed);
			}
		}
	}
fprintf (stderr, "----- MTF: %u    MTF SECOND LINE: %u\n",min_mtf, min_mtf_second_line);

fprintf (stderr, "HUFFMAN ON LITERALS: ");
print_clock(); start_clock();

	escaped = optimize_esc (graph, inbuf, in_len, min_no_esc, &min_start_esc, min_mtf, min_mtf_second_line);
        huffman_tree_len_bits = 0;
        huffed = generate_huffman_trees_and_codebook (graph, inbuf, in_len, huffman_tree, &huffman_tree_len_bits, table, head_array, 8-min_no_esc);

fprintf (stderr, "CALCULATE (AGAIN) USING FINALLY DETERMINED PARAMETERS: ");
print_clock(); start_clock();

	output_path (graph, inbuf, in_len, out_len, rle_occ, rle_rank, rank_len, min_start_esc, min_no_esc, min_mtf, min_mtf_second_line, min_lz_lsb, min_lz2_bits, table, head_array, huffman_tree, huffman_tree_len_bits );

fprintf(stderr, "-----== ENCODING STATS ==-----\n#ESC: 			%u\n#ESCAPED LIT:		%u\n#LZ OFFSET LSB:		%u\n#LZ2 OFFSET BITS	%u\n#RANKED RLEs		%u\nIN:			%u\nOUT: without header	%u\n     before huffman\n------------------------------\n",
	min_no_esc, escaped, min_lz_lsb, min_lz2_bits, rank_len, in_len, (min_no_bits + 7)>>3);

//	*out_len = (min_no_bits + 7) >> 3; // this ignores the header !!! delete !!!

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


uint16_t decode_huffman_tree (struct huffman_node *nodes, uint16_t *next_available_node, uint8_t no_symbol_bits, uint8_t leading_bits ) {

    // read bit
    if ( up_GetBits(1) ) {
	// a leaf at tree's end? read symbol's bits!
	uint16_t ret = leading_bits | up_GetBits(no_symbol_bits);
        nodes[ret].l = -1;
        nodes[ret].r = -1;
        // return leaf number
	return (ret);
    } else {
        uint16_t current_node = *next_available_node;
        *next_available_node = *next_available_node + 1;
        nodes[current_node].l = decode_huffman_tree (nodes, next_available_node, no_symbol_bits, leading_bits);
        nodes[current_node].r = decode_huffman_tree (nodes, next_available_node, no_symbol_bits, leading_bits);
        return (current_node);
    }
}

uint32_t build_huffman_trees (uint16_t no_esc, struct huffman_node *nodes, uint16_t heads[256]) {


	uint8_t no_trees = 0x01 << no_esc;

	uint8_t no_symbol_bits = 8 - no_esc;
	uint16_t head = 256;
	for (uint8_t i=0; i<no_trees; i++) {
		uint8_t leading_bits = i << no_symbol_bits;
		if ( up_GetBits(1) ) {
			// a tree to follow
			heads[i] = decode_huffman_tree (nodes, &head, no_symbol_bits, leading_bits);
		} else {
			// no tree
			heads[i] = 0xFFFF;
		}
	}
}

uint8_t huffman_256_decode ( struct huffman_node *nodes, uint16_t head) {

   /* decode the symbol */
        uint16_t n = head;
        while (nodes[n].l != -1) { /* more efficent maybe: n > 255 ? */
            if ( up_GetBits(1) )
                n = nodes[n].r;
            else
                n = nodes[n].l;
	}
        /* special case: no branches, only one leaf in the tree, just get one bit and do nothing about it */
        if (n == head) up_GetBits(1);
	return (n);
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
	// move to front
	uint8_t mtf = up_GetBits(1);
	uint8_t mtf_second_line = mtf?up_GetBits(8):0;

	// build huffman trees for LITerals
	uint16_t heads[256]; for (int i = 0; i<256; i++) heads[i]=0xFFFF;
	struct huffman_node nodes [511];
	build_huffman_trees ( no_esc, nodes, heads);

fprintf(stderr, "-----== DECODING STATS ==-----\n#ESC: 			%6u\n#LZ OFFSET LSB:		%6u\n#LZ2 OFFSET BITS	%6u\n#RANKED RLEs		%6u\nIN incl. header:	%6u\nOUT:			%6u\n------------------------------\n",
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
						// ESCaped (maybe mtf-encoded) LITeral
						// get the new ESCape code
						current_esc8 = up_GetBits(no_esc);
						current_esc = current_esc8 << (8 - no_esc);
						// fetch the rest of the mtf-encoded LITeral
						if (heads[first_fetch] != 0xFFFF) {
							first_fetch = (first_fetch<<(8-no_esc)) | huffman_256_decode (nodes, heads[first_fetch]);
						} else {
							first_fetch = (first_fetch<<(8-no_esc)) | up_GetBits(8-no_esc);
						}
						outbuf[cnt++] = mtf?move_to_front_decode (first_fetch, alphabet, mtf_second_line):first_fetch;
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
			if (heads[first_fetch] != 0xFFFF) {
				first_fetch = (first_fetch<<(8-no_esc)) | huffman_256_decode (nodes, heads[first_fetch]);
			} else {
				first_fetch = (first_fetch<<(8-no_esc)) | up_GetBits(8-no_esc);
			}
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
