#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <zlib.h>
#include <ctype.h>

#include "kseq.h"

#define PACKAGE_VERSION "2.0"

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })


/**bwa index struct**/
typedef uint64_t bwtint_t;
typedef unsigned char ubyte_t;
#define OCC_INTV_SHIFT 7
#define OCC_INTERVAL   (1LL<<OCC_INTV_SHIFT)
#define OCC_INTV_MASK  (OCC_INTERVAL - 1)

#define bwt_bwt(b, k) ((b)->bwt[((k)>>7<<4) + sizeof(bwtint_t) + (((k)&0x7f)>>4)])
#define bwt_occ_intv(b, k) ((b)->bwt + ((k)>>7<<4))

/* retrieve a character from the $-removed BWT string. Note that
 * bwt_t::bwt is not exactly the BWT string and therefore this macro is
 * called bwt_B0 instead of bwt_B */
#define bwt_B0(b, k) (bwt_bwt(b, k)>>((~(k)&0xf)<<1)&3)

typedef struct {
	bwtint_t primary; // S^{-1}(0), or the primary index of BWT
	bwtint_t L2[5]; // C(), cumulative count
	bwtint_t seq_len; // sequence length
	bwtint_t bwt_size; // size of bwt, about seq_len/4
	uint32_t *bwt; // BWT
	uint32_t cnt_table[256];
	int sa_intv;
	bwtint_t n_sa;
	bwtint_t *sa;
} bwt_t;

typedef struct {
	int64_t offset;
	int32_t len;
	int32_t n_ambs;
	uint32_t gi;
	char *name, *anno;
} bntann1_t;

typedef struct {
	int64_t offset;
	int32_t len;
	char amb;
} bntamb1_t;

typedef struct {
	int64_t l_pac;
	int32_t n_seqs;
	uint32_t seed;
	bntann1_t *anns; // n_seqs elements
	int32_t n_holes;
	bntamb1_t *ambs; // n_holes elements
	FILE *fp_pac;
} bntseq_t;

typedef struct {
	bwt_t    *bwt; // FM-index
	bntseq_t *bns; // information on the reference sequences
	uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
} bwaidx_t;


#define OUTPUT_ALL        0x01
#define OUTPUT_BEST       0x02

#define FA_FILE 0x01
#define FQ_FILE 0x02

#define SEED_SETP 0x00
#define EXON_STEP 0x01

#define KMER_FRONT 8
#define EXON_KMER 12

#define MAX_NAME_LENGTH 400
#define MAX_FILE_NUM 2000
#define MAX_STRING_LENGTH 2000000
#define MAX_READ_LENGTH 100000
#define MAX_CIGAR_LENGTH 200000

#define MAX_CHR_NUM 2000

struct file_name
{
    char name[MAX_NAME_LENGTH];
};
struct file_list
{
    struct file_name *file;
    int total;
};

struct chr_hash
{
    uint16_t back;
    unsigned int site;
};
struct chr_t
{
    char name[30];
    unsigned int length;
    unsigned int start_site;
    char *seq;
	int thread_num;
};
struct chr_list
{
    struct chr_t *list;
    int total;
};


struct splice_list
{
    unsigned int start,end;
    int num;
};
#define EXON_BUF_LENGTH 100000

#define READ_BUF_LENGTH 200//500//200//500//1000//500
#define SEED_BUF_LENGTH 500000
#define SEED_CAND_NUM 100
struct m_opt{
	char BWTpath[MAX_NAME_LENGTH];
	bwaidx_t *idx;
	char Hash_path[MAX_NAME_LENGTH];

	char Output_path[MAX_NAME_LENGTH];
	char gtf_file[MAX_NAME_LENGTH];

	struct file_list *input_file_1;

	unsigned int **c_hash;
	int *c_num;

	unsigned int **e_hash;
	int *e_num;

	int file_flag;
	struct chr_list *chr;

	int change_length;

	int area;
	int score_t;
	int match,miss,insert,gap,splice;

	int thread_num;
	int input_mode;

	int deep_mode;
	int pass;

	struct splice_list *SP;
	int SP_num;
	int SP_MAX;
	int SP_block_MAX;

	int64_t l_pac;

	int cut_tail;

	int thread_block;
	int result_block;

	int un_num;
};
struct m_opt *opt;
struct seed_t
{
    int start;
    int length;
    int chr_order;
    uint64_t pos;
    char flag;
    int score;
    uint64_t abs;

    int last;
    int lorder;
};
struct read_t
{
	char *name;
	int length;
	char *seq, *rseq, *qual,*rqual;
};
struct result_t
{
	int read_order;

	char *cigar;
	int chr;
	uint64_t ref_pos;
	int strand;
	int score;
};

struct snp_t
{
    unsigned int start;
    unsigned int end;
    int length;
    char type;
    uint16_t seq;
    int num;
    int dep;
};
struct snp_t *snp_buf;
struct snp_list_t
{
    struct snp_t *snp;
    uint64_t total;
};
struct job_seed
{
    FILE *in_file;
    FILE *out_file;

    FILE *exon_file;
    FILE *splice_file;

    FILE *un_file;
};
struct NW_list
{
    char *text;
    unsigned int *site;
};
struct DPC_t
{
    char cigar;
    int score;
    int m;
    int n;
};
struct max_t
{
    int score;
    int m;
};
struct Splice_DP_t
{
    struct DPC_t *cigarD;
    //struct DPC_t *cigarA;
    int** r;
    int** t;
    int** d;
    int** a;

    int** r1;
    int** t1;

    struct max_t *max_d;
    struct max_t *max_a;
};
#define MAX_DP_AREA 7000

struct seed_a
{
    bwtint_t start,end;
};
struct cigar_t
{
    int l;
    char c;
};
#define MAX_CIGAR_BUF 30000
#define MAX_MID_CIGAR_BUF 30

struct splice_t
{
    unsigned start;
    unsigned end;
    int length;

    struct cigar_t *cigar;
    int cigar_num;

    int score;
};

struct exon_seed_t
{
    int read_start,read_end;
    unsigned int ref_start,ref_end;
    int score;
    int flag;
    int front;

    struct cigar_t *cigar;
    int cigar_num;
    int c_start,c_end;
};
#define EXON_SEED 5000


extern char base2char[5];
extern unsigned char nst_nt4_table[256];
extern int usage();
//index.c
extern int index_main(int argc, char *argv[]);
extern int load_hash(struct m_opt *opt);
extern int hash_c_find(unsigned int **hash_front,int *hash_num,uint16_t base,unsigned int site);
extern uint16_t char2base_16(char *readseq,int length);

//bwa_index.c
extern void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol);
extern bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k);
extern int load_index(struct m_opt *opt);
extern void bwt_destroy(bwt_t *bwt);
extern void bns_destroy(bntseq_t *bns);
extern void bwa_idx_destroy(bwaidx_t *idx);

//seed.c
extern int seed_cmp(const void *a,const void *b);
extern int seed_pos_cmp(const void *a,const void *b);
extern int seed_cmp_r(const void *a,const void *b);
extern void insert_seed_start(struct seed_t *seed,int *seed_num,int init_num,int max,struct seed_t seed_r);
extern void find_seed_s(int step,int seed_length,struct read_t *read,struct seed_t *seed,int *seed_num,struct seed_a *seed_w);
extern void find_seed(struct read_t *read,struct seed_t *seed,int *seed_num,int *flag,struct seed_a *seed_w);
//generate_cigar.c
extern void out_read(struct read_t *read,struct result_t *result,int result_num,FILE *out);
extern void generate_cigar(struct read_t *read,struct seed_t *seed,int *seed_order,int seed_order_n,struct seed_t *seed_a,struct result_t *result,int *result_num,int read_order,FILE *out,FILE *un,struct NW_list *NW1,struct NW_list *NW2,struct NW_list *NW3,struct Splice_DP_t *DP,struct splice_list *SP,int *SP_num);
extern int SP_cmp(const void *a,const void *b);
extern int find_SP_Order(struct splice_list *SP,unsigned int start);
extern void Update_SP(struct splice_list *sp);
extern void Update_SP_s(struct splice_list *SP,int *SP_num,unsigned int MAX_N,struct splice_list *sp,int flag);
extern int splice_pos(int chr_order,uint64_t start,uint64_t end,char *seq,int length,struct cigar_t *cigar,int *cigar_num,int *score,int max,int mode,struct Splice_DP_t *DP);
extern int cigar_score(struct cigar_t *cigar,int cigar_num);
extern int Known_splice_pos(int chr_order,uint64_t start,uint64_t end,char *seq,int length,struct cigar_t *cigar,int *cigar_num,int *score,int max,struct Splice_DP_t *DP);

//DP.c
extern int Tail_DP(char *chr,unsigned int start,unsigned int end,int mode,char* string,int length,struct cigar_t *cigar,int *cigar_num,struct Splice_DP_t *DP,unsigned int *start_pos,int S);
extern int NW(char *ref,int n,char* text,int m,struct cigar_t *cigar,int *cigar_num,struct Splice_DP_t *DP);
extern int Area_DP(char *chr,unsigned int start,unsigned int end,char* string,int length,struct cigar_t *cigar,int *cigar_num,struct Splice_DP_t *DP);
extern int Splice_DP(char *chr,unsigned int start,unsigned int end,char* string,int length,struct cigar_t *cigar,int *cigar_num,struct Splice_DP_t *DP);
extern int tail_seed(int seed_length,int chr,char strand,uint64_t pos,char *seq,int length,int mode,struct seed_t *seed_o);

extern int seed_align(struct m_opt *opt);
//recheck.c
extern int re_check(struct m_opt *opt);


