//#include <stdlib.h>
//#include <ctype.h>
#include <stdio.h> // !!!
//#include <errno.h>
//#include <fcntl.h>
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>

struct code_t {
    uint8_t code[32];
    uint8_t code_len;
};

struct node {
    int32_t count, l, r;
    bool used;
    struct code_t code;
};



uint32_t huffman_tree_encode (struct node *nodes, uint32_t current_node, 
			      struct code_t current_code,
			      uint8_t *encoded_tree, uint32_t *tree_position ) {

    if (nodes[current_node].l == -1) {
	// write 1
	encoded_tree[*tree_position >> 3] = (encoded_tree[*tree_position >> 3] << 1) | 0x01;
	*tree_position = *tree_position + 1;
	uint8_t cn = current_node;
	// write 8-bit uncompressed value
	for (int i=7; i >= 0; i--) {
	    encoded_tree[*tree_position >> 3] = (encoded_tree[*tree_position >> 3] << 1) | ((cn >> i) & 0x01);
	    *tree_position = *tree_position + 1;
	}
	// complete missing shifts and write code to leaf
	for (int i=current_code.code_len; (i & 0x07) != 0; i++)
	    current_code.code[i >> 3] = current_code.code[i >> 3] << 1;
	for (int i=0; i < 32; i++)
	    nodes[current_node].code.code[i] = current_code.code[i];
	nodes[current_node].code.code_len = current_code.code_len;
    } else {
	// write (0lr) where l and r are the recursively determined respective sub-trees
	encoded_tree[*tree_position >> 3] = encoded_tree[*tree_position >> 3] << 1;
	*tree_position = *tree_position + 1;

	// construct codeword (l=0)
	current_code.code[current_code.code_len >> 3] = current_code.code[current_code.code_len >> 3] << 1;
	current_code.code_len++;
	// recursively parse left sub-tree
	huffman_tree_encode (nodes, nodes[current_node].l, current_code, encoded_tree, tree_position);

	// construct codeword (r=1)
	current_code.code_len--;
	current_code.code[current_code.code_len >> 3] = current_code.code[current_code.code_len >> 3] | 0x01;
	current_code.code_len++;
	// recursively parse right sub-tree
	huffman_tree_encode (nodes, nodes[current_node].r, current_code, encoded_tree, tree_position);
    }
}

uint32_t huffman_256_encode (uint8_t *outbuf,uint32_t *out_len,
                             const uint8_t *inbuf, uint16_t in_len) {

    if (in_len > 0) {
        struct node nodes[512] = {};

        /* count occurences */
        for (int i=0; i < in_len; i++) nodes[inbuf[i]].count++;

        /* mark the ones without hits as already used, so they
           they do not get checked later */
        for (int i=0; i < 256; i++) {
    	    if (nodes[i].count == 0) nodes[i].used = true;
	    nodes[i].l = -1;
	    nodes[i].r = -1;
        }

        uint32_t head = 255;
        do {
	    uint32_t min = in_len + 1;
	    uint32_t min_idx = -1;
	    head++;
            for (int i=0; i < head; i++) {
	        if (!nodes[i].used) {
	            if (nodes[i].count < min) {
		        min = nodes[i].count;
		        min_idx = i;
	            }
	        }
            }
            nodes[head].l = min_idx;
	    nodes[head].count = nodes[min_idx].count;
	    nodes[min_idx].used = true;

	    min = in_len + 1;
	    min_idx = -1;
            for (int i=0; i < head; i++) {
	        if (!nodes[i].used) {
	            if (nodes[i].count < min) {
		        min = nodes[i].count;
		        min_idx = i;
	            }
	        }
            }
	    if (min_idx != -1) {
                nodes[head].r = min_idx;
	        nodes[head].count += nodes[min_idx].count;
	        nodes[min_idx].used = true;
	    } else {
	        head=nodes[head].l;
	        break;
	    }
        } while (true);

        /* output tree and generate codebook */

        uint32_t max_out_len = *out_len << 3;
        *out_len=0;

        struct code_t code = {};

        huffman_tree_encode (nodes, head, code, outbuf, out_len );
/*    for (int i=*out_len; (i & 0x07) != 0; i++)
        outbuf[i >> 3] = outbuf[i >> 3] << 1;*/

/*        for (int i=0; i < *out_len; i++) {
            if (  ((outbuf[i >> 3] << (i & 0x07)) & 0x80) != 0 )
	    printf ("1"); else printf ("0");
        } +/

        printf ("\n SPACE: %u(%lx)\n",nodes[32].code.code_len, nodes[32].code.code[0]);
        printf ("\n NULL : %u(%lx)\n",nodes[0].code.code_len, nodes[0].code.code[0]);

        /* output data using code book */

//    uint32_t max_out_len =*out_len << 3;
//    *out_len = 0;

        // put data length
        for (int i=sizeof(in_len) * 8 -1; i >= 0; i--) {
	    outbuf[*out_len >> 3] = (outbuf[*out_len >> 3] << 1) | ((in_len >> i) & 0x01);
            *out_len = *out_len + 1;
        }

        for (int i=0; i < in_len; i++) {
	    /* special case: data consists of all the same byte value(s)
	       so the tree has no branches, just one leaf */
	    if (inbuf[i] == head) {
	        outbuf[*out_len >> 3] = outbuf[*out_len >> 3] << 1;
	        *out_len = *out_len + 1;
	    } else { /* regular case */
	        for (int j=0; j < nodes[inbuf[i]].code.code_len; j++) {
		    outbuf[*out_len >> 3] = (outbuf[*out_len >> 3] << 1) |
                                            ((nodes[inbuf[i]].code.code[j >> 3] >> (7 - (j & 0x07))) & 0x01);
		    *out_len = *out_len + 1;
	        }
	    }
        }
        for (int i=*out_len; (i & 0x07) != 0; i++)
            outbuf[i >> 3] = outbuf[i >> 3] << 1;

/*        for (int i=0; i < *out_len; i++) {
	if (  ((outbuf[i >> 3] << (i & 0x07)) & 0x80) != 0 )
	    printf ("1"); else printf ("0");
        }
        printf ("\n"); */

        *out_len = (*out_len+7) >> 3;
    } else {
	// empty inbuffer --> empty outbuffer
	*out_len = 0;
    }
}



uint32_t huffman_tree_decode (struct node *nodes, uint32_t *next_available_node,
			      uint8_t *encoded_tree, uint32_t *tree_position ) {

    // read bit
    if ( ((encoded_tree[*tree_position >> 3] << (*tree_position & 0x07)) & 0x80) != 0 ) {
	*tree_position = *tree_position + 1;
        // read 8 bits
	uint32_t ret = 0;
	for (int i=0; i < 8; i++) {
	    ret = (ret << 1) | (encoded_tree[*tree_position >> 3] >> (7 - (*tree_position & 0x07)) & 0x01);
	    *tree_position = *tree_position + 1;
	}
	// return leaf number
	nodes[ret].l = -1;
	nodes[ret].r = -1;
	return (ret);
    } else {
	*tree_position = *tree_position + 1;
	uint32_t current_node = *next_available_node;
	*next_available_node = *next_available_node + 1;

	nodes[current_node].l = huffman_tree_decode (nodes, next_available_node, encoded_tree, tree_position);
	nodes[current_node].r = huffman_tree_decode (nodes, next_available_node, encoded_tree, tree_position);
	return (current_node);
    }
}

uint32_t huffman_256_decode (uint8_t *outbuf,uint16_t *out_len,
                             uint8_t *inbuf, uint32_t in_len) {

    uint32_t max_out_len = *out_len;
    /* rebuild tree from input stream */
    struct node nodes[512] = {};
    uint32_t head = 256;
    uint32_t pos = 0;
    head = huffman_tree_decode (nodes, &head, inbuf, &pos);
    /* decode length */
    for (int i=sizeof(*out_len) * 8 -1; i >= 0; i--) {
	*out_len = (*out_len << 1) | ((inbuf[pos >> 3] >> (7 - (pos & 0x07))) & 0x01);
        pos = pos + 1;
    }
   /* decode the symbols */
   for (int i=0; i < *out_len; i++) {
        uint32_t n = head;
	while (nodes[n].l != -1) { /* more efficent maybe: n > 255 ? */
    	    if ( ((inbuf[pos >> 3] >> (7 - (pos & 0x07))) & 0x01) != 0)
	        n = nodes[n].r;
	    else
	        n = nodes[n].l;
            pos = pos + 1;
	}
	outbuf[i] = n;
	/* special case: no branches, only one leaf in the tree */
	if (n == head) pos = pos + 1;
    }
}



int main (void) {

    char input[] = "the best things in life are for free.";
//    char input[] = "";

    char output[100];
    uint32_t outlen = 100;

    huffman_256_encode (output, &outlen, input, sizeof(input));
    printf ("huffman: %u --> %u bytes\n",sizeof(input), outlen);

    char output2[100];
    uint16_t outlen2 = 100;

    huffman_256_decode (output2, &outlen2, output, outlen);
    return(0);
}

