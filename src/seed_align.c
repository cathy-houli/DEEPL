#include "main.h"

static pthread_mutex_t ReadLock,OutputLock;

int generate_flag_single(struct result_t *read)
{
    int flag = 0;

    if(read->cigar==NULL) {flag = 0x0004;return flag;}
    else
    {
        if(read->strand==1) flag|=0x0010;
    }
    return flag;
}
int read_file(struct read_t *read,FILE *IN_file,int mode)//0 FA 1 FQ
{
	int i = 0,j = 0;
    int l;
    char *f_line = (char *)calloc(MAX_STRING_LENGTH, sizeof(char));
    char *name = (char *)calloc(MAX_NAME_LENGTH, sizeof(char));

    if (opt->input_mode==FA_FILE)
    {
        while ((i<opt->thread_block)&&(fgets(f_line,MAX_STRING_LENGTH,IN_file)!=NULL))
        {
            sscanf(f_line,">%s",name);
            l = strlen(name);
            read[i].name = (char *)calloc(l+1,1);
            strncpy(read[i].name,name,l);
            read[i].name[l] = '\0';

            if (fgets(f_line,MAX_STRING_LENGTH,IN_file)!=NULL)
            {
                l = strlen(f_line)-1;
                read[i].length = l;
                read[i].seq = (char *)calloc(l+1,1);
                read[i].rseq = (char *)calloc(l+1,1);
                strncpy(read[i].seq,f_line,l);
                read[i].seq[l] = '\0';

                read[i].qual = NULL;
                read[i].rqual = NULL;

                for(j = 0;j<l;j++)
                {
                    read[i].seq[l-j-1] = toupper(read[i].seq[l-j-1]);
                    if((read[i].seq[l-j-1]=='A')||(read[i].seq[l-j-1]=='a')) read[i].rseq[j] = 'T';
                    else if((read[i].seq[l-j-1]=='C')||(read[i].seq[l-j-1]=='c')) read[i].rseq[j] = 'G';
                    else if((read[i].seq[l-j-1]=='G')||(read[i].seq[l-j-1]=='g')) read[i].rseq[j] = 'C';
                    else if((read[i].seq[l-j-1]=='T')||(read[i].seq[l-j-1]=='t')) read[i].rseq[j] = 'A';
                    else read[i].rseq[j] = 'N';
                }
                read[i].rseq[l] = '\0';
            }
            i++;
        }
    }
    else if (opt->input_mode==FQ_FILE)
    {
        while ((i<opt->thread_block)&&(fgets(f_line,MAX_STRING_LENGTH,IN_file)!=NULL))
        {
            sscanf(f_line,"@%s",name);
            l = strlen(name);
            read[i].name = (char *)calloc(l+1,1);
            strncpy(read[i].name,name,l);
            read[i].name[l] = '\0';

            if (fgets(f_line,MAX_STRING_LENGTH,IN_file)!=NULL)
            {
                l = strlen(f_line)-1;
                if(l>=MAX_READ_LENGTH)
                    printf("1\n");
                read[i].length = l;
                read[i].seq = (char *)calloc(l+1,1);
                read[i].rseq = (char *)calloc(l+1,1);
                strncpy(read[i].seq,f_line,l);
                read[i].seq[l] = '\0';

                for(j = 0;j<l;j++)
                {
                    read[i].seq[l-j-1] = toupper(read[i].seq[l-j-1]);
                    if((read[i].seq[l-j-1]=='A')||(read[i].seq[l-j-1]=='a')) read[i].rseq[j] = 'T';
                    else if((read[i].seq[l-j-1]=='C')||(read[i].seq[l-j-1]=='c')) read[i].rseq[j] = 'G';
                    else if((read[i].seq[l-j-1]=='G')||(read[i].seq[l-j-1]=='g')) read[i].rseq[j] = 'C';
                    else if((read[i].seq[l-j-1]=='T')||(read[i].seq[l-j-1]=='t')) read[i].rseq[j] = 'A';
                    else read[i].rseq[j] = 'N';
                }
                read[i].rseq[l] = '\0';
            }
            fgets(f_line,MAX_STRING_LENGTH,IN_file);
            if (fgets(f_line,MAX_STRING_LENGTH,IN_file)!=NULL)
            {
                read[i].qual= (char *)calloc(l+1,1);
                read[i].rqual = (char *)calloc(l+1,1);
                strncpy(read[i].qual,f_line,l);
                read[i].qual[l] = '\0';

                for(j = 0;j<l;j++)
                    read[i].rqual[j] = read[i].qual[l-j-1];
                read[i].rqual[l] = '\0';
            }
            i++;
        }
    }
    free(f_line);
    free(name);
    return i;
}
void free_read(struct read_t *read,int read_num)
{
    int i;
    for(i = 0;i<read_num;i++)
    {
        if (read->name!=NULL){free(read->name);read->name = NULL;}
        if (read->seq!=NULL) {free(read->seq);read->seq = NULL;}
        if (read->rseq!=NULL) {free(read->rseq);read->rseq = NULL;}
        if (read->qual!=NULL) {free(read->qual);read->qual = NULL;}
        if (read->rqual!=NULL) {free(read->rqual);read->rqual = NULL;}
    }
}
void out_read(struct read_t *read,struct result_t *result,int result_num,FILE *out)//,FILE *un)
{
    int i;
    int flag;
    for(i = 0;i<result_num;i++)
    {
        flag = generate_flag_single(&(result[i]));
        if(result[i].cigar==NULL)
        {
            //opt->un_num++;
            if (opt->input_mode==FA_FILE)
            fprintf(out,"%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t*\n",read[result[i].read_order].name,flag,read[result[i].read_order].seq);
            else if (opt->input_mode==FQ_FILE)
            fprintf(out,"%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",read[result[i].read_order].name,flag,read[result[i].read_order].seq,read[result[i].read_order].qual);
        }
        else
        {
            if (opt->input_mode==FA_FILE)
            fprintf(out,"%s\t%d\t%s\t%lu\t0\t%s\t*\t0\t0\t%s\t*\n",
                        read[result[i].read_order].name,flag,opt->chr->list[result[i].chr].name,result[i].ref_pos+1,result[i].cigar,
                        (result[i].strand==0)?(read[result[i].read_order].seq):(read[result[i].read_order].rseq));
            else if (opt->input_mode==FQ_FILE)
            fprintf(out,"%s\t%d\t%s\t%lu\t0\t%s\t*\t0\t0\t%s\t%s\n",
                        read[result[i].read_order].name,flag,opt->chr->list[result[i].chr].name,result[i].ref_pos+1,result[i].cigar,
                        (result[i].strand==0)?(read[result[i].read_order].seq):(read[result[i].read_order].rseq),(result[i].strand==0)?(read[result[i].read_order].qual):(read[result[i].read_order].rqual));
        }
        if(result[i].cigar!=NULL){free(result[i].cigar);result[i].cigar=NULL;}
    }
}

void extend_seed_tail_forward_double(int chr_order,char strand,char *seq,struct seed_t *seed,int *seed_num,uint64_t pos,int Sstart,int init_N)
{
    int i;
    struct seed_t seed_r[10];
    int flag = 0;

    if(Sstart>=4)
        flag = tail_seed(4,chr_order,strand,pos,seq,Sstart,0,seed_r);
    if(flag)
    {
        for(i = 0;i<flag;i++)
        {
            seed_r[i].abs = seed_r[i].pos-seed_r[i].start;
            seed_r[i].score = 0;
            seed_r[i].flag = 0;
            seed_r[i].last = chr_order;

            insert_seed_start(seed,seed_num,init_N,SEED_BUF_LENGTH,seed_r[i]);

            if(i==0)
                extend_seed_tail_forward_double(chr_order,strand,seq,seed,seed_num,seed_r[i].pos,seed_r[i].start,init_N);
        }
    }
}
void extend_seed_tail_backward_double(int chr_order,char strand,char *seq,int length,struct seed_t *seed,int *seed_num,uint64_t pos,int Sstart,int init_N)
{
    struct seed_t seed_r[10];
    int flag = 0;

    if(length-Sstart>=4)
        flag = tail_seed(4,chr_order,strand,pos-1,&seq[Sstart],length-Sstart,1,seed_r);
    if(flag)
    {
        int i;
        for(i = 0;i<flag;i++)
        {
            seed_r[i].start += Sstart;
            seed_r[i].abs = seed_r[i].pos-seed_r[i].start;
            seed_r[i].score = 0;
            seed_r[i].flag = 0;
            seed_r[i].last = chr_order;

            insert_seed_start(seed,seed_num,init_N,SEED_BUF_LENGTH,seed_r[i]);

            if(i==0)
                extend_seed_tail_backward_double(chr_order,strand,seq,length,seed,seed_num,seed_r[i].pos+seed_r[i].length,seed_r[i].start+seed_r[i].length,init_N);
        }
    }
}

struct best_t
{
    int seed_order;
    int score;
};
int best_cmp(const void *a,const void *b)
{
    struct best_t *EA,*EB;
    EA = (struct best_t *)a;
    EB = (struct best_t *)b;

    return EB->score-EA->score;
}
#define CAND_AREA 200000

int cut_tail(int *seed_order,int *seed_order_n,struct seed_t *seed)
{
    int i;
    uint64_t start;

    start = seed[seed_order[(*seed_order_n)-1]].pos;
    for(i = *seed_order_n-2;i>=0;i--)
    {
        if((seed[seed_order[i]].abs>seed[seed_order[i+1]].abs)&&(seed[seed_order[i]].abs-seed[seed_order[i+1]].abs>30))
        {
            if(seed[seed_order[i+1]].pos+seed[seed_order[i+1]].length-start<=opt->cut_tail)
                (*seed_order_n) = i+1;
            break;
        }
        else if(seed[seed_order[i]].pos+seed[seed_order[i]].length-start>opt->cut_tail) break;
    }

    start = seed[seed_order[0]].pos+seed[seed_order[0]].length;
    for(i = 0;i<(*seed_order_n)-1;i++)
    {
        if((seed[seed_order[i]].abs>seed[seed_order[i+1]].abs)&&(seed[seed_order[i]].abs-seed[seed_order[i+1]].abs>30))
        {
            if(start-seed[seed_order[i]].pos<=opt->cut_tail)
            {
                (*seed_order_n)-=i+1;
                return i+1;
            }
            break;
        }
        else if(start-seed[seed_order[i]].pos>opt->cut_tail) break;
    }

    return 0;
}
void find_cand(struct read_t *read,struct seed_t *seed,int *seed_num,struct seed_t *seed_a,struct best_t *best,struct result_t *result,int *result_num,int read_order,FILE *out,FILE *un,
               struct NW_list *NW1,struct NW_list *NW2,struct NW_list *NW3,struct Splice_DP_t *DP)
{
    int max_o = -1;
    int max_s = 0;
    int seed_order[10000];
    int seed_order_n = 0;

    int i = 0,j = 0;
    int x = 0;
    int chr_order = 0;
    unsigned int chr_start = 0;

    uint64_t LB,RB;
    //uint64_t LB;
    int init_num = *seed_num;
    char strand = 0;
    //extend
    for(i = 0; i<(*seed_num); i++)
    {
        if((seed[i].start == 0)&&(seed[i].pos == 0)&&(seed[i].length == 0)) continue;
        if(seed[i].pos<seed[i].start) {seed[i].start = 0;seed[i].pos = 0;seed[i].length = 0;seed[i].score = 0;seed[i].last = -1;continue;}

        if(seed[i].pos>=opt->l_pac)
        {
            strand = 1;
            for(j = 0;j<opt->chr->total;j++)
            {
                if((seed[i].pos >= (opt->l_pac<<1)-opt->chr->list[j].start_site-opt->chr->list[j].length)
                    &&(seed[i].pos < (opt->l_pac<<1)-opt->chr->list[j].start_site))
                    break;
            }
        }
        else{
            strand = 0;
            for(j = 0;j<opt->chr->total;j++)
            {
                if((seed[i].pos >= opt->chr->list[j].start_site)
                    &&(seed[i].pos <opt->chr->list[j].start_site+opt->chr->list[j].length))
                    break;
            }
        }
        if((strand == 0)&&(seed[i].pos+seed[i].length>=opt->chr->list[j].start_site+opt->chr->list[j].length))
        {seed[i].start = 0;seed[i].pos = 0;seed[i].length = 0;seed[i].score = 0;seed[i].last = -1;continue;}
        if((strand == 1)&&((opt->l_pac<<1)-seed[i].pos-seed[i].length<opt->chr->list[j].start_site))
        {seed[i].start = 0;seed[i].pos = 0;seed[i].length = 0;seed[i].score = 0;seed[i].last = -1;continue;}

        //seed_sxtend
        seed[i].last = j;
        chr_order = j;
        chr_start = opt->chr->list[j].start_site;
        int front = 0;
        front = seed[i].start;
        for(x = 0;x<front;x++)
        {
            if(seed[i].pos>=opt->l_pac)
            {
                if((seed[i].pos < (opt->l_pac<<1)-opt->chr->list[chr_order].start_site-opt->chr->list[chr_order].length)
                    ||(seed[i].pos >= (opt->l_pac<<1)-opt->chr->list[chr_order].start_site)) break;
                if ((3-nst_nt4_table[(int)opt->chr->list[chr_order].seq[(opt->l_pac<<1)-seed[i].pos-chr_start]])!=(nst_nt4_table[(int)read->seq[seed[i].start-1]])) break;
                else {seed[i].pos--;seed[i].start--;seed[i].length++;}
            }
            else
            {
                if((seed[i].pos < opt->chr->list[chr_order].start_site)
                    ||(seed[i].pos >=opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length))
                    break;
                if(nst_nt4_table[(int)opt->chr->list[chr_order].seq[seed[i].pos-1-chr_start]]!=nst_nt4_table[(int)read->seq[seed[i].start-1]])break;
                else {seed[i].pos--;seed[i].start--;seed[i].length++;}
            }
        }
        for(x = seed[i].start+seed[i].length;x<read->length;x++)
        {
            if(seed[i].pos>=opt->l_pac)
            {
                if((seed[i].pos+seed[i].length < (opt->l_pac<<1)-opt->chr->list[chr_order].start_site-opt->chr->list[chr_order].length)
                    ||(seed[i].pos+seed[i].length >= (opt->l_pac<<1)-opt->chr->list[chr_order].start_site)) break;
                if((3-nst_nt4_table[(int)opt->chr->list[chr_order].seq[(opt->l_pac<<1)-seed[i].pos-seed[i].length-1-chr_start]])!=(nst_nt4_table[(int)read->seq[x]])) break;
                else seed[i].length++;
            }
            else
            {
                if((seed[i].pos+seed[i].length < opt->chr->list[chr_order].start_site)
                    ||(seed[i].pos+seed[i].length >=opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length))
                    break;
                if(nst_nt4_table[(int)opt->chr->list[chr_order].seq[seed[i].pos+seed[i].length-chr_start]]!=nst_nt4_table[(int)read->seq[x]]) break;
                else seed[i].length++;
            }
        }
    }
    //unique and extend seed
    qsort(seed,*seed_num,sizeof(struct seed_t),seed_cmp);
    int start = 0;
    for(i = 0; i<init_num; i++)
    {
        while((seed[i].start >seed[start].start)&&(start<(*seed_num))) start++;
        for(j = start;j<i;j++)
        {
            if(seed[j].length == 0)continue;
            if((seed[i].start == seed[j].start)&&(seed[i].pos == seed[j].pos))
            {seed[i].start = 0;seed[i].pos = 0;seed[i].length = 0;seed[i].score = 0;seed[i].last = -1;break;}
        }
        if((seed[i].start == 0)&&(seed[i].pos == 0)&&(seed[i].length == 0)) continue;
    }
    qsort(seed,*seed_num,sizeof(struct seed_t),seed_cmp);
    for(i = 0; i<init_num; i++)
    {
        if(seed[i].length!=0)
        {
            if(seed[i].pos>=opt->l_pac) strand = 1;
            else strand = 0;

            chr_order = seed[i].last;
            extend_seed_tail_forward_double(chr_order,strand,read->seq,seed,seed_num,seed[i].pos,seed[i].start,init_num);
            extend_seed_tail_backward_double(chr_order,strand,read->seq,read->length,seed,seed_num,seed[i].pos+seed[i].length,seed[i].start+seed[i].length,init_num);
        }
    }

    qsort(seed,*seed_num,sizeof(struct seed_t),seed_cmp);
    start = 0;
    for(i = 0; i<(*seed_num); i++)
    {
        while((seed[i].start >seed[start].start)&&(start<(*seed_num))) start++;
        for(j = start;j<i;j++)
        {
            if(seed[j].length == 0)continue;
            if((seed[i].start == seed[j].start)&&(seed[i].pos == seed[j].pos))
            {seed[i].start = 0;seed[i].pos = 0;seed[i].length = 0;seed[i].score = 0;seed[i].last = -1;break;}
        }
    }

    //int contig_num = 0;
    //unique_seed(seed,*seed_num,read->length);
    qsort(seed,*seed_num,sizeof(struct seed_t),seed_pos_cmp);
    int k = 0;
    while((k<(*seed_num))&&(seed[k].start == 0)&&(seed[k].pos == 0)&&(seed[k].length == 0)) k++;
    start = k;
    for(i = 0; i<(*seed_num); i++)
    {
        if((seed[i].start == 0)&&(seed[i].pos == 0)&&(seed[i].length == 0)) continue;

        seed[i].abs = seed[i].pos-seed[i].start;
        RB=seed[i].pos+seed[i].length;//opt->change_length;//Insert
        //LB
        if(seed[i].pos>=opt->l_pac)
            LB = max(((opt->l_pac<<1)-opt->chr->list[seed[i].last].start_site-opt->chr->list[seed[i].last].length),(seed[i].pos>CAND_AREA?seed[i].pos-CAND_AREA:0));
        else
            LB = max(opt->chr->list[seed[i].last].start_site,(seed[i].pos>CAND_AREA?seed[i].pos-CAND_AREA:0));

        //从头开始对每个seed赋分，即本seed长度+前置最优路径长度
        max_o = -1;
        max_s = -1000;
        seed[i].last = -1;
        seed[i].lorder = 0;
        seed[i].score = 0;
        while((LB >seed[start].pos)&&(start<(*seed_num))) start++;
        for(j = start;j<i;j++)
        {
            if(seed[j].length == 0)continue;
            //if((seed[i].start == seed[j].start)&&(seed[i].pos == seed[j].pos))
            //{seed[i].start = 0;seed[i].pos = 0;seed[i].length = 0;seed[i].score = 0;seed[i].last = -1;break;}

            if((seed[i].start > seed[j].start)
                &&((seed[i].start+seed[i].length)>(seed[j].start+seed[j].length))
                &&(
                  ((seed[i].abs > seed[j].abs)?((seed[i].abs-seed[j].abs)<seed[j].length):((seed[j].abs-seed[i].abs)<seed[i].length))
                  ||((RB>=(seed[j].pos+seed[i].start-seed[j].start))&&(LB<=(seed[j].pos+seed[i].start-seed[j].start)))

                  )

                //&&(RB>=(seed[j].pos+seed[i].start-seed[j].start))
                //&&(LB<=(seed[j].pos+seed[i].start-seed[j].start))
             )
            {
                if(seed[i].start<=seed[j].start+seed[j].length)
                    max_s = seed[j].score + (seed[i].start+seed[i].length-seed[j].start-seed[j].length)*opt->match;
                else
                    max_s = seed[j].score + seed[i].length*opt->match;

                if(seed[i].abs<seed[j].abs)
                    max_s -=(seed[j].abs-seed[i].abs)*opt->match+opt->gap;//max_s -=(seed[j].abs-seed[i].abs)*opt->match+4+(seed[j].abs-seed[i].abs-1)*opt->gap;//max_s -=(seed[j].abs-seed[i].abs)*opt->match+opt->gap; //max_s -=(seed[j].abs-seed[i].abs)*opt->match+4+(seed[j].abs-seed[i].abs-1)*opt->gap;//max_s -=(seed[j].abs-seed[i].abs)+2;
                else if((seed[i].abs>seed[j].abs)&&(seed[i].abs-seed[j].abs<=opt->change_length))
                    max_s -=opt->gap;//max_s -=4+(seed[i].abs-seed[j].abs-1)*opt->gap;//max_s -=opt->gap;//max_s -=4+(seed[i].abs-seed[j].abs-1)*opt->gap;//max_s -=opt->gap;
                else if((seed[i].abs>seed[j].abs)&&(seed[i].abs-seed[j].abs>opt->change_length))
                    max_s -=opt->splice;

                if(seed[i].score <= max_s)
                {
                    max_o = j;
                    seed[i].score = max_s;
                }
            }
        }
        if(max_o == -1)
            seed[i].score = seed[i].length*opt->match;
        else
            seed[i].last = max_o;
    }
    //回溯
    max_o = 0;
    max_s = -1000;
    int flag = 0;
    int max_e = 0;
    init_num = (*result_num);

    for(i = (*seed_num)-1; i>=0; i--)
    {
        best[i].seed_order = i;
        best[i].score = seed[i].score;
    }
    qsort(best,(*seed_num),sizeof(struct best_t),best_cmp);

    if(best[0].score>=(float)0.3*read->length*opt->match)
    {
        max_o = best[0].seed_order;
        seed_order_n = 0;
        while(max_o!=-1)
        {
            //if(seed[max_o].length>=8)
            {
                seed_order[seed_order_n] = max_o;
                seed_order_n++;
            }
            seed[max_o].flag = 1;
            if(seed_order_n>=600) {flag = 1;break;}

            if(seed[max_o].last!=-1)
                max_o = seed[max_o].last;
            else max_o = -1;
        }
        if(flag == 0)
        {
            max_s=cut_tail(seed_order,&seed_order_n,seed);
            generate_cigar(read,seed,&(seed_order[max_s]),seed_order_n,seed_a,result,result_num,read_order,out,un,NW1,NW2,NW3,DP);
        }
    }
    else if((opt->deep_mode==1)&&(best[0].score>=(float)0.2*read->length*opt->match))//else if((opt->deep_mode==1)&&(opt->pass==2)&&(best[0].score>=(float)0.2*read->length*opt->match))
    {
        max_o = best[0].seed_order;
        seed_order_n = 0;
        while(max_o!=-1)
        {
            //if(seed[max_o].length>=8)
            {
                seed_order[seed_order_n] = max_o;
                seed_order_n++;
            }
            seed[max_o].flag = 1;
            if(seed_order_n>=600) {flag = 1;break;}

            if(seed[max_o].last!=-1)
                max_o = seed[max_o].last;
            else max_o = -1;
        }
        if(flag == 0)
        {
            max_s=cut_tail(seed_order,&seed_order_n,seed);
            generate_cigar(read,seed,&(seed_order[max_s]),seed_order_n,seed_a,result,result_num,read_order,out,un,NW1,NW2,NW3,DP);
        }
    }
    else
    {
        result[(*result_num)].read_order = read_order;
        result[(*result_num)].cigar = NULL;
        result[(*result_num)].chr = 0;
        result[(*result_num)].ref_pos  = 0;
        result[(*result_num)].strand = 0;
        (*result_num)++;
        if((*result_num)>=opt->result_block)
        {
            pthread_mutex_lock(&OutputLock);
            out_read((read-read_order),result,*result_num,out);
            pthread_mutex_unlock(&OutputLock);
            (*result_num) = 0;
        }
    }

    for(i = 1; i<(*seed_num); i++)
    {
        if(best[i].score>=(float)0.7*read->length*opt->match)
        {
            max_e = best[i].seed_order;
            seed_order_n = 0;
            flag = 0;
            while(max_e!=-1)
            {
                if(seed[max_e].flag==1) {flag = 1;break;}

                //if(seed[max_e].length>=8)
                {
                seed_order[seed_order_n] = max_e;
                seed_order_n++;
                }
                seed[max_e].flag = 1;

                if(seed_order_n>=600) {flag = 1;break;}

                if(seed[max_e].last!=-1)
                    max_e = seed[max_e].last;
                else max_e = -1;
            }
            if(flag == 0)
            {
                max_s=cut_tail(seed_order,&seed_order_n,seed);
                generate_cigar(read,seed,&(seed_order[max_s]),seed_order_n,seed_a,result,result_num,read_order,out,un,NW1,NW2,NW3,DP);
            }
        }
        else break;
    }
    if(init_num == (*result_num))
    {
        result[(*result_num)].read_order = read_order;
        result[(*result_num)].cigar = NULL;
        result[(*result_num)].chr = 0;
        result[(*result_num)].ref_pos  = 0;
        result[(*result_num)].strand = 0;
        (*result_num)++;
        if((*result_num)>=opt->result_block)
        {
            pthread_mutex_lock(&OutputLock);
            out_read((read-read_order),result,*result_num,out);
            pthread_mutex_unlock(&OutputLock);
            (*result_num) = 0;
        }
    }
}
void *seed_align_core(void* arg)
{
    struct job_seed *job = (struct job_seed *)arg;

	struct seed_t *seed = (struct seed_t *)malloc((SEED_BUF_LENGTH+1)*sizeof(struct seed_t));
	int seed_num;
	struct seed_t *seed_a = (struct seed_t *)malloc((SEED_BUF_LENGTH+1)*sizeof(struct seed_t));
	struct best_t *best = (struct best_t *)malloc((SEED_BUF_LENGTH+1)*sizeof(struct best_t));

	struct read_t *read = (struct read_t *)malloc(opt->thread_block*sizeof(struct read_t));
	int read_num = 0;

	struct result_t *result = (struct result_t *)malloc((opt->result_block+1)*sizeof(struct result_t));
	int result_num = 0;

	int *flag = (int *)calloc(MAX_READ_LENGTH+5, sizeof(int));
    struct seed_a *seed_w = (struct seed_a *)calloc(MAX_READ_LENGTH+5, sizeof(struct seed_a));

	struct NW_list *NW1 = (struct NW_list *)calloc(1,sizeof(struct NW_list));
	NW1->text = (char *)malloc((MAX_READ_LENGTH+1)*sizeof(char));
    NW1->site = (unsigned int *)malloc((MAX_READ_LENGTH+1)* sizeof(unsigned int));

    struct NW_list *NW2 = (struct NW_list *)calloc(1,sizeof(struct NW_list));
	NW2->text = (char *)malloc((MAX_READ_LENGTH+1)*sizeof(char));
    NW2->site = (unsigned int *)malloc((MAX_READ_LENGTH+1)* sizeof(unsigned int));

    struct NW_list *NW3 = (struct NW_list *)calloc(1,sizeof(struct NW_list));
	NW3->text = (char *)malloc((MAX_READ_LENGTH+1)*sizeof(char));
    NW3->site = (unsigned int *)malloc((MAX_READ_LENGTH+1)* sizeof(unsigned int));

    struct Splice_DP_t *DP = (struct Splice_DP_t *)calloc(1,sizeof(struct Splice_DP_t));
    DP->cigarD = (struct DPC_t *)malloc((2*MAX_DP_AREA)*sizeof(struct DPC_t));
    //DP->cigarA = (struct DPC_t *)malloc((2*MAX_DP_AREA)*sizeof(struct DPC_t));

    DP-> r = (int **)malloc((MAX_DP_AREA+1)*sizeof(int *));
	DP-> t = (int **)malloc((MAX_DP_AREA+1)*sizeof(int *));
	DP-> r1 = (int **)malloc((MAX_DP_AREA+1)*sizeof(int *));
	DP-> t1 = (int **)malloc((MAX_DP_AREA+1)*sizeof(int *));
	DP-> d = (int **)malloc((MAX_DP_AREA+1)*sizeof(int *));
	DP-> a = (int **)malloc((MAX_DP_AREA+1)*sizeof(int *));
	DP-> max_d = (struct max_t *)malloc((MAX_DP_AREA+1)*sizeof(struct max_t));
	DP-> max_a = (struct max_t *)malloc((MAX_DP_AREA+1)*sizeof(struct max_t));

	int i;
	for (i = 0; i <=MAX_DP_AREA; i++)
	{
		DP->r[i] = (int *)malloc((MAX_DP_AREA+1)*sizeof(int));
		DP->t[i] = (int *)malloc((MAX_DP_AREA+1)*sizeof(int));
		DP->r1[i] = (int *)malloc((MAX_DP_AREA+1)*sizeof(int));
		DP->t1[i] = (int *)malloc((MAX_DP_AREA+1)*sizeof(int));
		DP->d[i] = (int *)malloc((MAX_DP_AREA+1)*sizeof(int));
		DP->a[i] = (int *)malloc((MAX_DP_AREA+1)*sizeof(int));
	}

	int break_flag = 0;

	while (1)
	{
        //read_file;
        pthread_mutex_lock(&ReadLock);
        read_num = read_file(read,job->in_file,opt->input_mode);
		while (read_num == 0)
		{
            if(opt->file_flag >= (opt->input_file_1->total-1)){break_flag = 1;break;}
            else
            {
                opt->file_flag++;
                fclose(job->in_file);
                job->in_file = fopen(opt->input_file_1->file[opt->file_flag].name,"r");
                read_num = read_file(read,job->in_file,opt->input_mode);
            }
		}
		pthread_mutex_unlock(&ReadLock);

		if(break_flag) break;
        for (i = 0;i<read_num;i++)
        {
            seed_num = 0;
            find_seed(read+i,seed,&seed_num,flag,seed_w);
            if(opt->deep_mode==1)
                find_seed_s(5,15,read+i,seed,&seed_num,seed_w);
                //find_seed_r(read+i,seed,&seed_num,flag,seed_w);
            else
                find_seed_s(25,25,read+i,seed,&seed_num,seed_w);

            if(seed_num!=0)
                find_cand(read+i,seed,&seed_num,seed_a,best,result,&result_num,i,job->out_file,job->un_file,NW1,NW2,NW3,DP);
        }

        pthread_mutex_lock(&OutputLock);
        out_read(read,result,result_num,job->out_file);
        result_num = 0;
        pthread_mutex_unlock(&OutputLock);

        free_read(read,read_num);
	}

	free(seed);
	free(seed_a);
	free(best);
	free(read);
	free(result);

	free(NW1->text);
    free(NW1->site);
    free(NW1);

	free(NW2->text);
    free(NW2->site);
    free(NW2);

	free(NW3->text);
    free(NW3->site);
    free(NW3);

    for (i = 0; i <=MAX_DP_AREA; i++)
	{
		free(DP->r[i]); free(DP->t[i]);free(DP->r1[i]); free(DP->t1[i]); free(DP->d[i]); free(DP->a[i]);
	}
	free(DP->r); free(DP->t);free(DP->r1); free(DP->t1); free(DP->d); free(DP->a);
	free(DP->cigarD);free(DP-> max_d);free(DP-> max_a);
	free(DP);

    free(flag);
    free(seed_w);

    return (void*)(0);
}
int seed_align(struct m_opt *opt)
{
    struct job_seed *job = (struct job_seed *)calloc(1,sizeof(struct job_seed));
    int i;

    job->in_file = fopen(opt->input_file_1->file[0].name,"r");
    if (job->in_file == NULL){return 1;}

    job->out_file = fopen(opt->Output_path,"w");
    if (job->out_file == NULL){return 1;}

	pthread_t *pthreads = malloc(sizeof(pthread_t) * opt->thread_num);

    for (i = 0; i < opt->thread_num; i++) pthread_create(&pthreads[i], NULL, seed_align_core, job);
	for (i = 0; i < opt->thread_num; i++) pthread_join(pthreads[i], NULL);

    fclose(job->in_file);
    fclose(job->out_file);

    free(pthreads);
    free(job);

    return 0;
}

