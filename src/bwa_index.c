#include "main.h"

bntseq_t *bns_restore_core(const char *ann_filename, const char* amb_filename, const char* pac_filename)
{
	char str[1024];
	FILE *fp;
	const char *fname;
	bntseq_t *bns;
	long long xx;
	int i;
	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	{ // read .ann
		fp = fopen(fname = ann_filename, "r");
		fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);
		bns->l_pac = xx;
		bns->anns = (bntann1_t*)calloc(bns->n_seqs, sizeof(bntann1_t));
		for (i = 0; i < bns->n_seqs; ++i) {
			bntann1_t *p = bns->anns + i;
			char *q = str;
			int c;
			// read gi and sequence name
			fscanf(fp, "%u%s", &p->gi, str);
			p->name = strdup(str);
			// read fasta comments
			while (str - q < (int)(sizeof(str) - 1) && (c = fgetc(fp)) != '\n' && c != EOF) *q++ = c;
			while (c != '\n' && c != EOF) c = fgetc(fp);

			*q = 0;
			if (q - str > 1) p->anno = strdup(str + 1); // skip leading space
			else p->anno = strdup("");
			// read the rest
			fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
			p->offset = xx;
		}
		fclose(fp);
	}
	{ // read .amb
		//int64_t l_pac;
		int32_t n_seqs;
		fp = fopen(fname = amb_filename, "r");
		fscanf(fp, "%lld%d%d", &xx, &n_seqs, &bns->n_holes);
		//l_pac = xx;
		//xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
		bns->ambs = bns->n_holes? (bntamb1_t*)calloc(bns->n_holes, sizeof(bntamb1_t)) : 0;
		for (i = 0; i < bns->n_holes; ++i) {
			bntamb1_t *p = bns->ambs + i;
			fscanf(fp, "%lld%d%s", &xx, &p->len, str);
			p->offset = xx;
			p->amb = str[0];
		}
		fclose(fp);
	}
	{ // open .pac
		bns->fp_pac = fopen(pac_filename, "rb");
	}
	return bns;
}
bntseq_t *bns_restore(const char *prefix)
{
	char ann_filename[256], amb_filename[256], pac_filename[256];
	strcat(strcpy(ann_filename, prefix), ".ann");
	strcat(strcpy(amb_filename, prefix), ".amb");
	strcat(strcpy(pac_filename, prefix), ".pac");
	return bns_restore_core(ann_filename, amb_filename, pac_filename);
}

static bwtint_t fread_fix(FILE *fp, bwtint_t size, void *a)
{
	const unsigned int bufsize = 0x1000000; // 16M block
	bwtint_t offset = 0;
	while (size) {
		bwtint_t x = bufsize < size ? bufsize : size;
		x = fread(((char*)a + offset), 1, x, fp);
		size -= x; offset += x;
	}
	return offset;
}
void bwt_restore_sa(const char *fn, bwt_t *bwt)
{
	char skipped[256];
	FILE *fp;
	bwtint_t primary;

	fp = fopen(fn, "rb");
	fread(&primary, sizeof(bwtint_t), 1, fp);
	//xassert(primary == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
	fread(skipped, sizeof(bwtint_t), 4, fp); // skip
	fread(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	fread(&primary, sizeof(bwtint_t), 1, fp);
	//xassert(primary == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");

	bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
	bwt->sa[0] = -1;

	fread_fix(fp, sizeof(bwtint_t) * (bwt->n_sa - 1), bwt->sa + 1);
	fclose(fp);
}
void bwt_gen_cnt_table(bwt_t *bwt)
{
	int i, j;
	for (i = 0; i != 256; ++i) {
		uint32_t x = 0;
		for (j = 0; j != 4; ++j)
			x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
		bwt->cnt_table[i] = x;
	}
}
bwt_t *bwt_restore_bwt(const char *fn)
{
	bwt_t *bwt;
	FILE *fp;

	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	fp = fopen(fn, "rb");
	fseek(fp, 0, SEEK_END);
	bwt->bwt_size = (ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
	bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4);
	fseek(fp, 0, SEEK_SET);
	fread(&bwt->primary, sizeof(bwtint_t), 1, fp);
	fread(bwt->L2 + 1, sizeof(bwtint_t), 4, fp);
	fread_fix(fp, bwt->bwt_size << 2, bwt->bwt);
	bwt->seq_len = bwt->L2[4];
	fclose(fp);
	bwt_gen_cnt_table(bwt);

	return bwt;
}
bwt_t *bwa_idx_load_bwt(const char *hint)
{
	char *tmp;
	bwt_t *bwt;

	tmp = (char*)calloc(strlen(hint) + 5, 1);
	strcat(strcpy(tmp, hint), ".bwt"); // FM-index
	bwt = bwt_restore_bwt(tmp);
	strcat(strcpy(tmp, hint), ".sa");  // partial suffix array (SA)
	bwt_restore_sa(tmp, bwt);
	free(tmp);

	return bwt;
}
void bwt_destroy(bwt_t *bwt)
{
	if (bwt == 0) return;
	free(bwt->sa); free(bwt->bwt);
	free(bwt);
}
void bns_destroy(bntseq_t *bns)
{
	if (bns == 0) return;
	else {
		int i;
		if (bns->fp_pac) fclose(bns->fp_pac);
		free(bns->ambs);
		for (i = 0; i < bns->n_seqs; ++i) {
			free(bns->anns[i].name);
			free(bns->anns[i].anno);
		}
		free(bns->anns);
		free(bns);
	}
}
static inline int __occ_aux(uint64_t y, int c)
{
	// reduce nucleotide counting to bits counting
	y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
	// count the number of 1s in y
	y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}
bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c)
{
	bwtint_t n;
	uint32_t *p, *end;

	if (k == bwt->seq_len) return bwt->L2[c+1] - bwt->L2[c];
	if (k == (bwtint_t)(-1)) return 0;
	k -= (k >= bwt->primary); // because $ is not in bwt

	// retrieve Occ at k/OCC_INTERVAL
	n = ((bwtint_t*)(p = bwt_occ_intv(bwt, k)))[c];
	p += sizeof(bwtint_t); // jump to the start of the first BWT cell

	// calculate Occ up to the last k/32
	end = p + (((k>>5) - ((k&~OCC_INTV_MASK)>>5))<<1);
	for (; p < end; p += 2) n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);

	// calculate Occ
	n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
	if (c == 0) n -= ~k&31; // corrected for the masked bits

	return n;
}
// an analogy to bwt_occ() but more efficient, requiring k <= l
void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol)
{
	bwtint_t _k, _l;
	_k = (k >= bwt->primary)? k-1 : k;
	_l = (l >= bwt->primary)? l-1 : l;
	if (_l/OCC_INTERVAL != _k/OCC_INTERVAL || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
		*ok = bwt_occ(bwt, k, c);
		*ol = bwt_occ(bwt, l, c);
	} else {
		bwtint_t m, n, i, j;
		uint32_t *p;
		if (k >= bwt->primary) --k;
		if (l >= bwt->primary) --l;
		n = ((bwtint_t*)(p = bwt_occ_intv(bwt, k)))[c];
		p += sizeof(bwtint_t);
		// calculate *ok
		j = k >> 5 << 5;
		for (i = k/OCC_INTERVAL*OCC_INTERVAL; i < j; i += 32, p += 2)
			n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
		m = n;
		n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
		if (c == 0) n -= ~k&31; // corrected for the masked bits
		*ok = n;
		// calculate *ol
		j = l >> 5 << 5;
		for (; i < j; i += 32, p += 2)
			m += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
		m += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~l&31)<<1)) - 1), c);
		if (c == 0) m -= ~l&31; // corrected for the masked bits
		*ol = m;
	}
}
static inline bwtint_t bwt_invPsi(const bwt_t *bwt, bwtint_t k) // compute inverse CSA
{
	bwtint_t x = k - (k > bwt->primary);
	x = bwt_B0(bwt, x);
	x = bwt->L2[x] + bwt_occ(bwt, k, x);
	return k == bwt->primary? 0 : x;
}
bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k)
{
	bwtint_t sa = 0, mask = bwt->sa_intv - 1;
	while (k & mask) {
		++sa;
		k = bwt_invPsi(bwt, k);
	}
	/* without setting bwt->sa[0] = -1, the following line should be
	   changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
	return sa + bwt->sa[k/bwt->sa_intv];
}
struct job_ref
{
    struct m_opt *opt;
    int i;
    unsigned int length;
};
void *read_reference(void* arg)
{
    int i = 0;
    unsigned int j = 0,site = 0;
    int base;
    struct job_ref *job = (struct job_ref *)arg;
    int start = -1,end = -1;

    for(i = 0;i<job->opt->chr->total;i++)
    {
        if (job->opt->chr->list[i].start_site >= (job->i * job->length))
        {if(start==-1) start = i;}

        if (job->opt->chr->list[i].start_site >= ((job->i + 1) * job->length)) break;
        end = i;
    }
	for (i = start ; i <=end; i++)
	{
        job->opt->chr->list[i].seq = (char *)calloc(job->opt->chr->list[i].length, sizeof(char));
        for(j = 0;j<job->opt->chr->list[i].length;j++)
        {
            site = job->opt->chr->list[i].start_site+j;
            base = opt->idx->pac[site >> 2] >> ((~site & 3) << 1) & 3;
            switch (base)
            {
            case 0: job->opt->chr->list[i].seq[j] = 'A'; break;
            case 1: job->opt->chr->list[i].seq[j] = 'C'; break;
            case 2: job->opt->chr->list[i].seq[j] = 'G'; break;
            case 3: job->opt->chr->list[i].seq[j] = 'T'; break;
            default:job->opt->chr->list[i].seq[j] = 'N';
            }
        }
	}
	return (void*)(0);
}
void bwa_idx_load(struct m_opt *opt,const char *hint)
{
	fprintf(stderr, "Load the genome index files...");
	opt->idx->bwt = bwa_idx_load_bwt(hint);//load bwt
	opt->idx->bns = bns_restore(hint);// load reference informatiom

	opt->idx->pac = (uint8_t*)calloc(opt->idx->bns->l_pac/4+1, 1);// load reference pac seq
	fseek(opt->idx->bns->fp_pac, 0, SEEK_SET);
	fread(opt->idx->pac, 1, opt->idx->bns->l_pac / 4 + 1, opt->idx->bns->fp_pac);

	int i = 0;
	unsigned int start = 0;
	for (i = 0; i < opt->idx->bns->n_seqs; i++) //turn pac to char seq
	{
        opt->chr->list[i].length = opt->idx->bns->anns[i].len;
        strcpy(opt->chr->list[i].name,opt->idx->bns->anns[i].name);
        opt->chr->list[i].start_site = start;
        opt->chr->list[i].thread_num = 0;
        start += opt->idx->bns->anns[i].len;
	}
	opt->chr->total = opt->idx->bns->n_seqs;

	struct job_ref *job = (struct job_ref *)calloc(opt->thread_num,sizeof(struct job_ref));
	pthread_t *pthreads = malloc(sizeof(pthread_t) * opt->thread_num);
    for (i = 0; i < opt->thread_num; i++)
	{
        job[i].opt = opt;
		job[i].i = i;
		job[i].length = start/opt->thread_num;
		pthread_create(&pthreads[i], NULL, read_reference, job + i);
	}
	for (i = 0; i < opt->thread_num; i++) pthread_join(pthreads[i], NULL);

	free(job);
	free(pthreads);

	fprintf(stderr, "\n");
}
void bwa_idx_destroy(bwaidx_t *idx)
{
	if (idx->bwt) bwt_destroy(idx->bwt);
	if (idx->bns) bns_destroy(idx->bns);
	if (idx->pac) free(idx->pac);
}
int check_bwt_index(char *index)
{
	char file[MAX_NAME_LENGTH];

	sprintf(file,"%s.ann",index);
	if((access(file,0))== -1) return 1;

	sprintf(file,"%s.amb",index);
	if((access(file,0))== -1) return 1;

	sprintf(file,"%s.pac",index);
	if((access(file,0))== -1) return 1;

	return 0;
}
int load_index(struct m_opt *opt)
{
    if ((opt->BWTpath[0] != '\0') && !check_bwt_index(opt->BWTpath)) bwa_idx_load(opt,opt->BWTpath);
    else
    {
        fprintf(stderr, "Error! Cannot get reference index...\n");
        return 1;
    }
    return 0;
}
