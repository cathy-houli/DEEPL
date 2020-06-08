#include "main.h"

int seed_cmp(const void *a,const void *b)
{
    struct seed_t *EA,*EB;
    EA = (struct seed_t *)a;
    EB = (struct seed_t *)b;

    if(EA->start==EB->start)
    {
        if (EA->pos==EB->pos) return 0;
        else if (EA->pos<EB->pos) return -1;
        else  return 1;
    }
    else
    {
        if(EA->start < EB->start) return -1;
        else return 1;
    }
}
int seed_pos_cmp(const void *a,const void *b)
{
    struct seed_t *EA,*EB;
    EA = (struct seed_t *)a;
    EB = (struct seed_t *)b;

    if(EA->pos==EB->pos)
    {
        if (EA->start==EB->start) return 0;
        else if (EA->start<EB->start) return -1;
        else  return 1;
    }
    else
    {
        if(EA->pos < EB->pos) return -1;
        else return 1;
    }
}
int seed_cmp_r(const void *a,const void *b)
{
    struct seed_t *EA,*EB;
    EA = (struct seed_t *)a;
    EB = (struct seed_t *)b;

    if(EA->start==EB->start)
    {
        if (EA->pos==EB->pos) return 0;
        else if (EA->pos>EB->pos) return -1;
        else  return 1;
    }
    else
    {
        if(EA->start < EB->start) return -1;
        else return 1;
    }
}
int find_seed_order_p(struct seed_t *seed,int seed_num,uint64_t pos)
{
    int left = 0, right = seed_num-1, middle;

    if (right == -1)
        return 0;
    while (left <= right)
    {
        middle = (left + right)/2;

        if (seed[middle].pos == pos)
        {
            while ((middle-1>=0)&&(seed[middle-1].pos == pos))
                middle--;
            return middle;
        }
        else if (seed[middle].pos > pos)
            right = middle -1;
        else
            left = middle + 1;
    }
    return left;
}
void insert_seed_pos(struct seed_t *seed,int *seed_num,int init_num,int max,struct seed_t seed_r)
{
    int i = 0;
    int flag = 1;
    int order = find_seed_order_p(seed,init_num,seed_r.pos);
    for(i = order;i<*seed_num;i++)
    {
        if((seed[i].abs==seed_r.abs)&&(seed_r.pos==seed[i].pos)&&(seed_r.start=seed[i].start))
        {
            if((seed_r.length)>(seed[i].length))
                seed[i].length=seed_r.length;
            flag = 0;
            break;
        }
        else if(seed_r.pos<seed[i].pos) break;
    }
    if((flag)&&((*seed_num)<max))
        {memcpy(seed+(*seed_num),&seed_r,sizeof(struct seed_t));(*seed_num)++;}
}
int find_seed_order_s(struct seed_t *seed,int seed_num,int start)
{
    int left = 0, right = seed_num-1, middle;

    if (right == -1)
        return 0;
    while (left <= right)
    {
        middle = (left + right)/2;

        if (seed[middle].start == start)
        {
            while ((middle-1>=0)&&(seed[middle-1].start == start))
                middle--;
            return middle;
        }
        else if (seed[middle].start > start)
            right = middle -1;
        else
            left = middle + 1;
    }
    return left;
}
void insert_seed_start(struct seed_t *seed,int *seed_num,int init_num,int max,struct seed_t seed_r)
{
    int i = 0;
    int flag = 1;
    int order = find_seed_order_s(seed,init_num,seed_r.start);
    for(i = order;i<*seed_num;i++)
    {
        if((seed[i].abs==seed_r.abs)&&(seed_r.pos==seed[i].pos)&&(seed_r.start=seed[i].start))
        {
            if((seed_r.length)>(seed[i].length))
                seed[i].length=seed_r.length;
            flag = 0;
            break;
        }
        else if(seed_r.pos<seed[i].pos) break;
    }
    if((flag)&&((*seed_num)<max))
        {memcpy(seed+(*seed_num),&seed_r,sizeof(struct seed_t));(*seed_num)++;}
}
void insert_seed1(struct seed_t *seed,int *seed_num,int init_num,int max,struct seed_t seed_r)
{
    int i = 0;
    int flag = 1;
    for(i = 0;i<*seed_num;i++)
    {
        if((seed[i].start == 0)&&(seed[i].pos == 0)&&(seed[i].length == 0)) continue;
        if((seed[i].abs==seed_r.abs)&&(seed_r.start>=seed[i].start)&&(seed_r.start<=seed[i].start+seed[i].length))
        {
            if((seed_r.start+seed_r.length)>(seed[i].start+seed[i].length))
                seed[i].length+=seed_r.start+seed_r.length-seed[i].start-seed[i].length;
            flag = 0;
            break;
        }
    }
    if((flag)&&((*seed_num)<max))
        {memcpy(seed+(*seed_num),&seed_r,sizeof(struct seed_t));(*seed_num)++;}
}
void find_seed(struct read_t *read,struct seed_t *seed,int *seed_num,int *flag,struct seed_a *seed_w)
{
    int i = 0,j = 0,x = 0;//back_step = 0;;

    int c = 0;
    bwtint_t k, l, ok, ol;

    for (i = 0;i<read->length;i++)
    {
        //if(i%50==0) flag[i]=1;
        //else
        flag[i] = 0;
    }
    flag[read->length-1] = 1;
    //flag[read->length-3] = 1;
    //flag[read->length-2] = 1;
    //flag[read->length-4] = 1;
    flag[read->length-5] = 1;
    //flag[read->length-9] = 1;


    for (i = read->length-1;i>15;i--)
    {
        if(flag[i])
        {
            //bwt_find_seed;
            k = 0; l = opt->idx->bwt->seq_len;
            for (j = i;j>=0;j--)
            {
                if(nst_nt4_table[(int)read->seq[j]]>=4) {flag[i-1]=1;break;}
                c = nst_nt4_table[(int)read->seq[j]];
                bwt_2occ(opt->idx->bwt, k - 1, l, c, &ok, &ol);
                k = opt->idx->bwt->L2[c] + ok + 1;
                l = opt->idx->bwt->L2[c] + ol;
                seed_w[j].start = k;
                seed_w[j].end = l;

                if (k > l)
                {
                    if(((seed_w[j+1].end-seed_w[j+1].start+1) >SEED_CAND_NUM)||(i-j<16)) {flag[i-1]=1;break;}

                    for (x = 0;x<=seed_w[j+1].end-seed_w[j+1].start;x++)
                    {
                        seed[*seed_num].start = j+1;
                        seed[*seed_num].length = i-j;
                        seed[*seed_num].pos = bwt_sa(opt->idx->bwt,seed_w[j+1].start + x);
                        seed[*seed_num].last = -1;
                        seed[*seed_num].flag = 0;
                        seed[*seed_num].score = 0;
                        if((*seed_num)<SEED_BUF_LENGTH) (*seed_num)++;
                    }
                    if(j>=16)
                    {
                        flag[j]=1;
                        //flag[j-1]=1;
                        //flag[j-2]=1;
                        //flag[j-3]=1;
                        flag[j-4]=1;
                        //flag[j-8]=1;
                    }
                    //flag[j+1]=1;
                    //flag[j+2]=1;
                    //flag[j+3]=1;
                    flag[j+4]=1;
                    //flag[j+8]=1;
                    break;
                }
            }
            if(j == -1)
            {
                if (k <= l)
                {
                    if(((seed_w[j+1].end-seed_w[j+1].start+1) >SEED_CAND_NUM)||(i-j<16)) {flag[i-1]=1;continue;}

                    for (x = 0;x<=seed_w[j+1].end-seed_w[j+1].start;x++)
                    {
                        seed[*seed_num].start = j+1;
                        seed[*seed_num].length = i-j;
                        seed[*seed_num].pos = bwt_sa(opt->idx->bwt,seed_w[j+1].start + x);
                        seed[*seed_num].last = -1;
                        seed[*seed_num].flag = 0;
                        seed[*seed_num].score = 0;
                        if((*seed_num)<SEED_BUF_LENGTH)(*seed_num)++;
                    }
                }
            }
        }
    }
    //qsort(seed,*seed_num,sizeof(struct seed_t),seed_cmp);
}
void find_seed_r(struct read_t *read,struct seed_t *seed,int *seed_num,int *flag,struct seed_a *seed_w)
{
    int i = 0,j = 0,x = 0;//back_step = 0;;

    int c = 0;
    bwtint_t k, l, ok, ol;

    for (i = 0;i<read->length;i++)
    {
        //if(i%50==0) flag[i]=1;
        //else
            flag[i] = 0;
    }
    flag[read->length-1] = 1;
    //flag[read->length-3] = 1;
    //flag[read->length-2] = 1;
    //flag[read->length-4] = 1;
    flag[read->length-5] = 1;
    //flag[read->length-9] = 1;

    for (i = read->length-1;i>15;i--)
    {
        if(flag[i])
        {
            //bwt_find_seed;
            k = 0; l = opt->idx->bwt->seq_len;
            for (j = i;j>=0;j--)
            {
                if(nst_nt4_table[(int)read->rseq[j]]>=4) {flag[i-1]=1;break;}
                c = nst_nt4_table[(int)read->rseq[j]];
                bwt_2occ(opt->idx->bwt, k - 1, l, c, &ok, &ol);
                k = opt->idx->bwt->L2[c] + ok + 1;
                l = opt->idx->bwt->L2[c] + ol;
                seed_w[j].start = k;
                seed_w[j].end = l;

                if (k > l)
                {
                    if(((seed_w[j+1].end-seed_w[j+1].start+1) >SEED_CAND_NUM)||(i-j<16)) {flag[i-1]=1;break;}

                    for (x = 0;x<=seed_w[j+1].end-seed_w[j+1].start;x++)
                    {
                        seed[*seed_num].start = read->length-1-i;
                        seed[*seed_num].length = i-j;
                        seed[*seed_num].pos = (opt->l_pac<<1)-seed[*seed_num].length-bwt_sa(opt->idx->bwt,seed_w[j+1].start + x);
                        seed[*seed_num].last= -1;
                        seed[*seed_num].flag = 0;
                        seed[*seed_num].score = 0;
                        if((*seed_num)<SEED_BUF_LENGTH) (*seed_num)++;
                    }
                    if(j>=16)
                    {
                        flag[j]=1;
                        //flag[j-1]=1;
                        //flag[j-2]=1;
                        //flag[j-3]=1;
                        flag[j-4]=1;
                        //flag[j-8]=1;
                    }
                    //flag[j+1]=1;
                    //flag[j+2]=1;
                    //flag[j+3]=1;
                    flag[j+4]=1;
                    //flag[j+8]=1;
                    break;
                }
            }
            if(j == -1)
            {
                if (k <= l)
                {
                    if(((seed_w[j+1].end-seed_w[j+1].start+1) >SEED_CAND_NUM)||(i-j<16)) {flag[i-1]=1;continue;}

                    for (x = 0;x<=seed_w[j+1].end-seed_w[j+1].start;x++)
                    {
                        seed[*seed_num].start = read->length-1-i;
                        seed[*seed_num].length = i-j;
                        seed[*seed_num].pos = (opt->l_pac<<1)-seed[*seed_num].length-bwt_sa(opt->idx->bwt,seed_w[j+1].start + x);
                        seed[*seed_num].last = -1;
                        seed[*seed_num].flag = 0;
                        seed[*seed_num].score = 0;
                        if((*seed_num)<SEED_BUF_LENGTH) (*seed_num)++;
                    }
                }
            }
        }
    }
    //qsort(seed,*seed_num,sizeof(struct seed_t),seed_cmp);
}
void find_seed_s(int step,int seed_length,struct read_t *read,struct seed_t *seed,int *seed_num,struct seed_a *seed_w)
{
    int i = 0,j = 0,x = 0;//back_step = 0;;

    int c = 0;
    bwtint_t k, l, ok, ol;

    //struct seed_t seed_a;

    for (i = i*step;i*step+seed_length<read->length;i++)
    {
            //bwt_find_seed;
            k = 0; l = opt->idx->bwt->seq_len;
            for (j =i*step;j<i*step+seed_length;j++)
            {
                if(nst_nt4_table[(int)read->seq[j]]>=4) break;
                c = 3-nst_nt4_table[(int)read->seq[j]];
                bwt_2occ(opt->idx->bwt, k - 1, l, c, &ok, &ol);
                k = opt->idx->bwt->L2[c] + ok + 1;
                l = opt->idx->bwt->L2[c] + ol;
                seed_w[j].start = k;
                seed_w[j].end = l;

                if (k > l)
                    break;
            }
            if(j == i*step+seed_length)
            {
                if (k <= l)
                {
                    if((seed_w[j-1].end-seed_w[j-1].start+1) >SEED_CAND_NUM) continue;

                    for (x = 0;x<=seed_w[j-1].end-seed_w[j-1].start;x++)
                    {
                        seed[*seed_num].start = i*step;
                        seed[*seed_num].length = seed_length;
                        seed[*seed_num].pos = (opt->l_pac<<1)-seed[*seed_num].length-bwt_sa(opt->idx->bwt,seed_w[j-1].start + x);
                        seed[*seed_num].last = -1;
                        seed[*seed_num].score = 0;
                        seed[*seed_num].flag = 0;
                        seed[*seed_num].abs = seed[*seed_num].pos-seed[*seed_num].start;
                        if((*seed_num)<SEED_BUF_LENGTH) (*seed_num)++;
                        //insert_seed(seed,seed_num,SEED_BUF_LENGTH,seed_a);
                    }
                }
            }
    }
    //qsort(seed,*seed_num,sizeof(struct seed_t),seed_cmp);
}

