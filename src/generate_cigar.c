#include "main.h"
static pthread_mutex_t OutputLock,UpdateLock;

int SP_cmp(const void *a,const void *b)
{
    struct splice_list *EA,*EB;
    EA = (struct splice_list *)a;
    EB = (struct splice_list *)b;

    if(EA->start==EB->start)
    {
        if (EA->end==EB->end)
        {
            return 0;
        }
        else if (EA->end < EB->end) return -1;
        else  return 1;
    }
    else if(EA->start < EB->start) return -1;
    else return 1;
}
int find_SP_Order(struct splice_list *SP,unsigned int start)
{
    int left = 0, right = opt->SP_num-1, middle;

    if (right == -1)
        return -1;
    while (left <= right)
    {
        middle = (left + right)/2;

        if(SP[middle].start==start)
        {
            while((SP[middle-1].start==start)&&(middle>0)) middle--;
            return middle;
        }
        else if (SP[middle].start>start)
        {
            right = middle -1;
        }
        else
            left = middle + 1;
    }
    return left;
}
void Update_SP(struct splice_list *sp)
{
    //pthread_mutex_lock(&UpdateLock);
    if(opt->SP_num+1>=opt->SP_MAX)
    {
        opt->SP_MAX+=10000;
        opt->SP = (struct splice_list *)realloc(opt->SP,opt->SP_MAX*sizeof(struct splice_list));
    }
    memcpy(opt->SP+opt->SP_num,sp,1*sizeof(struct splice_list));
    opt->SP_num ++;
    //pthread_mutex_unlock(&UpdateLock);
}
void free_SP(int SN,struct splice_t *splice)
{
    int i = 0;
    for(i = 0;i<SN;i++)
    {
        if(splice[i].cigar!=NULL)free(splice[i].cigar);
    }
}
int splice_pos(int chr_order,uint64_t start,uint64_t end,char *seq,int length,struct cigar_t *cigar,int *cigar_num,int *score,int max,int mode,struct Splice_DP_t *DP)//MODE:1 GT-AG 2 CT-AC
{
    unsigned int GT[200];
    unsigned int AG[200];
    unsigned int CT[200];
    unsigned int AC[200];

    int GTN = 0,AGN = 0,CTN = 0,ACN = 0;

    unsigned int i = 0,j;
    unsigned int Sstart,Send;

    int TH = 5;

    Sstart = start+length+TH;
    Send = end-length-TH;

    char *chr = opt->chr->list[chr_order].seq;

    char *text = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
    char *ref = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
    for(j = 0;j<length;j++)
        text[j] = seq[j];
    text[length] = '\0';

    struct cigar_t *t_cigar = (struct cigar_t *)calloc(MAX_CIGAR_BUF,sizeof(struct cigar_t));
    int t_num = 0;

    struct cigar_t *b_cigar = (struct cigar_t *)calloc(MAX_CIGAR_BUF,sizeof(struct cigar_t));

    int istart;
    int length_r;

    struct splice_t best_splice;
    best_splice.score = -10000;
    best_splice.cigar = b_cigar;
    best_splice.cigar_num = 0;

    struct splice_t temp_splice;
    temp_splice.cigar = b_cigar;
    int SN = 0;
    int k = 0;

    if(mode==0)
    {
    for(i = start;i<Sstart;i++)
    {
        if((nst_nt4_table[(int)chr[i]]==2)&&(nst_nt4_table[(int)chr[i+1]]==3))
        {GT[GTN] = i;GTN++;if(GTN>=200){(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(b_cigar);free(t_cigar);return 0;}}
    }
    for(i = Send;i<=end+1;i++)
    {
        if((nst_nt4_table[(int)chr[i-2]]==0)&&(nst_nt4_table[(int)chr[i-1]]==2))
        {AG[AGN] = i;AGN++;if(AGN>=200){(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(b_cigar);free(t_cigar);return 0;}}
    }
    for(k = 0;k<GTN;k++)
    {
        for(j = 0;j<AGN;j++)
        {
            temp_splice.start = GT[k];
            temp_splice.end = AG[j];
            temp_splice.length = temp_splice.end-temp_splice.start;
            if(temp_splice.end-opt->change_length<=temp_splice.start) continue;

            length_r = (temp_splice.start-start)+(end+1-temp_splice.end);
            if(length_r>=MAX_READ_LENGTH) continue;
            for(i = start;i<temp_splice.start;i++)
                ref[i-start] = chr[i];
            istart = i-start;
            for(i = temp_splice.end;i<end+1;i++)
                ref[istart+i-temp_splice.end] = chr[i];
            ref[length_r] = '\0';

            SN++;
            if(SN>200) break;

            if(length_r==0)
            {
                temp_splice.cigar_num = 1;
                t_cigar[0].c ='I';
                t_cigar[0].l =length;
                t_num = 1;
                temp_splice.score = -4-(t_cigar[0].l-1)*opt->gap;
            }
            else
            {
                t_num = 0;
                temp_splice.score = NW(ref,length_r,text,length,t_cigar,&t_num,DP);
                temp_splice.cigar_num = t_num;
            }

            if((SN==1)||(temp_splice.score>best_splice.score))
            {
                memcpy(&best_splice,&temp_splice,1*sizeof(struct splice_t));
                memcpy(b_cigar,t_cigar,t_num*sizeof(struct cigar_t));
            }
        }
    }
    }
    else if(mode==1)
    {
    for(i = start;i<Sstart;i++)
    {
        if((nst_nt4_table[(int)chr[i]]==1)&&(nst_nt4_table[(int)chr[i+1]]==3))
        {CT[CTN] = i;CTN++;if(CTN>=200){(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(b_cigar);free(t_cigar);return 0;}}
    }
    for(i = Send;i<=end+1;i++)
    {
        if((nst_nt4_table[(int)chr[i-2]]==0)&&(nst_nt4_table[(int)chr[i-1]]==1))
        {AC[ACN] = i;ACN++;if(ACN>=200){(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(b_cigar);free(t_cigar);return 0;}}
    }
    for(k = 0;k<CTN;k++)
    {
        for(j = 0;j<ACN;j++)
        {
            temp_splice.start = CT[k];
            temp_splice.end = AC[j];
            temp_splice.length = temp_splice.end-temp_splice.start;
            if(temp_splice.end-opt->change_length<=temp_splice.start) continue;

            length_r = (temp_splice.start-start)+(end+1-temp_splice.end);
            if(length_r>=MAX_READ_LENGTH) continue;
            for(i = start;i<temp_splice.start;i++)
                ref[i-start] = chr[i];
            istart = i-start;
            for(i = temp_splice.end;i<end+1;i++)
                ref[istart+i-temp_splice.end] = chr[i];
            ref[length_r] = '\0';

            SN++;
            if(SN>200) break;

            if(length_r==0)
            {
                temp_splice.cigar_num = 1;
                t_cigar[0].c ='I';
                t_cigar[0].l =length;
                t_num = 1;
                temp_splice.score = -4-(t_cigar[0].l-1)*opt->gap;
            }
            else
            {
                t_num = 0;
                temp_splice.score = NW(ref,length_r,text,length,t_cigar,&t_num,DP);
                temp_splice.cigar_num = t_num;
            }

            if((SN==1)||(temp_splice.score>best_splice.score))
            {
                memcpy(&best_splice,&temp_splice,1*sizeof(struct splice_t));
                memcpy(b_cigar,t_cigar,t_num*sizeof(struct cigar_t));
            }
        }
    }
    }

    if((SN!=0)&&(best_splice.score!=-10000))
    {
        (*score) = best_splice.score;
        uint64_t ref_pos = start;
        int read_pos = 0;

        for(i = 0;i<best_splice.cigar_num;i++)
        {
            if((b_cigar[i].c=='M')||(b_cigar[i].c=='X')||(b_cigar[i].c=='D'))
            {
                if(ref_pos+b_cigar[i].l<=best_splice.start)
                {
                    cigar[(*cigar_num)].c = b_cigar[i].c;
                    cigar[(*cigar_num)].l = b_cigar[i].l;
                    (*cigar_num)++;
                    if((*cigar_num)>=max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(b_cigar);free(t_cigar);return 0;}

                    ref_pos+=b_cigar[i].l;
                    read_pos+=b_cigar[i].l;
                }
                else
                {
                    if(best_splice.start-ref_pos>0)
                    {
                        cigar[(*cigar_num)].c = b_cigar[i].c;
                        cigar[(*cigar_num)].l = best_splice.start-ref_pos;
                        (*cigar_num)++;
                        if((*cigar_num)>=max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(b_cigar);free(t_cigar);return 0;}
                        b_cigar[i].l-= best_splice.start-ref_pos;
                    }
                    break;
                }
            }
            else if(b_cigar[i].c=='I')
            {
                cigar[(*cigar_num)].c = b_cigar[i].c;
                cigar[(*cigar_num)].l = b_cigar[i].l;
                (*cigar_num)++;
                if((*cigar_num)>=max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(b_cigar);free(t_cigar);return 0;}
            }
            else {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(b_cigar);free(t_cigar);return 0;}
        }
        cigar[(*cigar_num)].c = 'N';
        cigar[(*cigar_num)].l = best_splice.length;
        (*cigar_num)++;
        if((*cigar_num)>=max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(b_cigar);free(t_cigar);return 0;}
        for(j = i;j<best_splice.cigar_num;j++)
        {
            cigar[(*cigar_num)].c = b_cigar[j].c;
            cigar[(*cigar_num)].l = b_cigar[j].l;
            (*cigar_num)++;
            if((*cigar_num)>=max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(b_cigar);free(t_cigar);return 0;}
        }
    }
    else {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free(b_cigar);return 0;}

    free(ref);
    free(text);
    free(t_cigar);
    free(b_cigar);

    return 0;
}
int Known_splice_pos(int chr_order,uint64_t start,uint64_t end,char *seq,int length,struct cigar_t *cigar,int *cigar_num,int *score,int max,struct Splice_DP_t *DP)//MODE:1 GT-AG 2 CT-AC
{
    unsigned int i = 0,j;
    unsigned int Sstart,Send;

    int TH = 5;

    Sstart = start+length+TH;
    Send = end-length-TH;

    char *chr = opt->chr->list[chr_order].seq;
    unsigned int  chr_start= opt->chr->list[chr_order].start_site;

    char *text = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
    char *ref = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
    for(j = 0;j<length;j++)
        text[j] = seq[j];
    text[length] = '\0';

    struct cigar_t *t_cigar = (struct cigar_t *)calloc(MAX_CIGAR_BUF,sizeof(struct cigar_t));
    int t_num = 0;

    struct cigar_t *b_cigar = (struct cigar_t *)calloc(MAX_CIGAR_BUF,sizeof(struct cigar_t));

    int istart;
    int length_r;

    struct splice_t best_splice;
    best_splice.score = -10000;
    best_splice.cigar = b_cigar;
    best_splice.cigar_num = 0;

    struct splice_t temp_splice;
    temp_splice.cigar = b_cigar;
    int SN = 0;
    int k = 0;

    int order = find_SP_Order(opt->SP,start+chr_start);
    if(order==-1) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free(b_cigar);return 0;}

    for(k = order;k<opt->SP_num;k++)
    {
        if(opt->SP[k].start>Sstart+chr_start) break;
        if(opt->SP[k].end<Send+chr_start) continue;
        if(opt->SP[k].end>end+chr_start) continue;

        {
            temp_splice.start = opt->SP[k].start-chr_start;
            temp_splice.end = opt->SP[k].end-chr_start;
            temp_splice.length = temp_splice.end-temp_splice.start;

            length_r = (temp_splice.start-start)+(end+1-temp_splice.end);
            if(length_r>=MAX_READ_LENGTH) continue;
            ref[length_r] = '\0';
            for(i = start;i<temp_splice.start;i++)
                ref[i-start] = chr[i];
            istart = i-start;
            for(i = temp_splice.end;i<end+1;i++)
                ref[istart+i-temp_splice.end] = chr[i];

            SN++;
            if(SN>200) break;

            if(length_r==0)
            {
                temp_splice.cigar_num = 1;
                t_cigar[0].c ='I';
                t_cigar[0].l =length;
                t_num = 1;
                temp_splice.score = -4-(t_cigar[0].l-1)*opt->gap;
            }
            else
            {
                t_num = 0;
                temp_splice.score = NW(ref,length_r,text,length,t_cigar,&(t_num),DP);
                temp_splice.score+=4+opt->SP[k].num;
                temp_splice.cigar_num = t_num;
            }
            if((SN==1)||(temp_splice.score>best_splice.score))
            {
                memcpy(&best_splice,&temp_splice,1*sizeof(struct splice_t));
                memcpy(b_cigar,t_cigar,t_num*sizeof(struct cigar_t));
            }
        }
    }

    if((SN!=0)&&(best_splice.score!=-10000))
    {
        (*score) = best_splice.score;
        uint64_t ref_pos = start;
        int read_pos = 0;

        for(i = 0;i<best_splice.cigar_num;i++)
        {
            if((b_cigar[i].c=='M')||(b_cigar[i].c=='X')||(b_cigar[i].c=='D'))
            {
                if(ref_pos+b_cigar[i].l<=best_splice.start)
                {
                    cigar[(*cigar_num)].c = b_cigar[i].c;
                    cigar[(*cigar_num)].l = b_cigar[i].l;
                    (*cigar_num)++;
                    if((*cigar_num)>max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free(b_cigar);return 0;}

                    ref_pos+=b_cigar[i].l;
                    read_pos+=b_cigar[i].l;
                }
                else
                {
                    if(best_splice.start-ref_pos>0)
                    {
                        cigar[(*cigar_num)].c = b_cigar[i].c;
                        cigar[(*cigar_num)].l = best_splice.start-ref_pos;
                        (*cigar_num)++;
                        if((*cigar_num)>max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free(b_cigar);return 0;}
                        b_cigar[i].l-= best_splice.start-ref_pos;
                    }
                    break;
                }
            }
            else if(b_cigar[i].c=='I')
            {
                cigar[(*cigar_num)].c = b_cigar[i].c;
                cigar[(*cigar_num)].l = b_cigar[i].l;
                (*cigar_num)++;
                if((*cigar_num)>max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free(b_cigar);return 0;}
            }
            else {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free(b_cigar);return 0;}
        }
        cigar[(*cigar_num)].c = 'N';
        cigar[(*cigar_num)].l = best_splice.length;
        (*cigar_num)++;
        if((*cigar_num)>max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free(b_cigar);return 0;}
        for(j = i;j<best_splice.cigar_num;j++)
        {
            cigar[(*cigar_num)].c = b_cigar[j].c;
            cigar[(*cigar_num)].l = b_cigar[j].l;
            (*cigar_num)++;
            if((*cigar_num)>max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free(b_cigar);return 0;}
        }
    }
    else {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free(b_cigar);return 0;}
    free(ref);
    free(text);
    free(t_cigar);
    free(b_cigar);
    return 0;
}
int splice_pos1(int chr_order,uint64_t start,uint64_t end,char *seq,int length,struct cigar_t *cigar,int *cigar_num,int *score,int max,int mode,struct Splice_DP_t *DP)//MODE:1 GT-AG 2 CT-AC
{
    unsigned int GT[200];
    unsigned int AG[200];
    unsigned int CT[200];
    unsigned int AC[200];

    int GTN = 0,AGN = 0,CTN = 0,ACN = 0;

    unsigned int i = 0,j;
    unsigned int Sstart,Send;

    int TH = 5;

    Sstart = start+length+TH;
    Send = end-length-TH;

    char *chr = opt->chr->list[chr_order].seq;

    char *text = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
    char *ref = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
    for(j = 0;j<length;j++)
        text[j] = seq[j];
    text[length] = '\0';

    struct cigar_t *t_cigar = (struct cigar_t *)calloc(MAX_CIGAR_BUF,sizeof(struct cigar_t));
    int t_num = 0;

    int istart;
    int length_r;

    struct splice_t splice[200];
    int SN = 0;
    int k = 0;

    if(mode==0)
    {
    for(i = start;i<Sstart;i++)
    {
        if((nst_nt4_table[(int)chr[i]]==2)&&(nst_nt4_table[(int)chr[i+1]]==3))
        {GT[GTN] = i;GTN++;if(GTN>=200){(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free_SP(SN,splice);free(t_cigar);return 0;}}
    }
    for(i = Send;i<=end+1;i++)
    {
        if((nst_nt4_table[(int)chr[i-2]]==0)&&(nst_nt4_table[(int)chr[i-1]]==2))
        {AG[AGN] = i;AGN++;if(AGN>=200){(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free_SP(SN,splice);free(t_cigar);return 0;}}
    }
    for(k = 0;k<GTN;k++)
    {
        for(j = 0;j<AGN;j++)
        {
            splice[SN].start = GT[k];
            splice[SN].end = AG[j];
            splice[SN].length = splice[SN].end-splice[SN].start;
            //splice[SN].site_score = splice_socre(0,1,chr,GT[k]-5)+splice_socre(0,2,chr,AG[j]-5);
            if(splice[SN].end-opt->change_length<=splice[SN].start) continue;
            if(SN>=200)continue;

            length_r = (splice[SN].start-start)+(end+1-splice[SN].end);
            if(length_r>=MAX_READ_LENGTH) continue;
            for(i = start;i<splice[SN].start;i++)
                ref[i-start] = chr[i];
            istart = i-start;
            for(i = splice[SN].end;i<end+1;i++)
                ref[istart+i-splice[SN].end] = chr[i];
            ref[length_r] = '\0';

            if(length_r==0)
            {
                splice[SN].cigar_num = 1;
                splice[SN].cigar = (struct cigar_t *)malloc(1*sizeof(struct cigar_t));
                splice[SN].cigar[0].c ='I';
                splice[SN].cigar[0].l =length;
                splice[SN].score = -4-(splice[SN].cigar[0].l-1)*opt->gap;
                SN++;
                if(SN>=200){(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free_SP(SN,splice);free(t_cigar);return 0;}
            }
            else
            {
                t_num = 0;
                splice[SN].score = NW(ref,length_r,text,length,t_cigar,&t_num,DP);
                splice[SN].cigar = (struct cigar_t *)malloc(t_num*sizeof(struct cigar_t));
                memcpy(splice[SN].cigar,t_cigar,t_num*sizeof(struct cigar_t));
                splice[SN].cigar_num = t_num;
                SN++;
                if(SN>=200){(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free_SP(SN,splice);free(t_cigar);return 0;}
            }

        }
    }
    }
    else if(mode==1)
    {
    for(i = start;i<Sstart;i++)
    {
        if((nst_nt4_table[(int)chr[i]]==1)&&(nst_nt4_table[(int)chr[i+1]]==3))
        {CT[CTN] = i;CTN++;if(CTN>=200){(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free_SP(SN,splice);free(t_cigar);return 0;}}
    }
    for(i = Send;i<=end+1;i++)
    {
        if((nst_nt4_table[(int)chr[i-2]]==0)&&(nst_nt4_table[(int)chr[i-1]]==1))
        {AC[ACN] = i;ACN++;if(ACN>=200){(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free_SP(SN,splice);free(t_cigar);return 0;}}
    }
    for(k = 0;k<CTN;k++)
    {
        for(j = 0;j<ACN;j++)
        {
            splice[SN].start = CT[k];
            splice[SN].end = AC[j];
            splice[SN].length = splice[SN].end-splice[SN].start;
            //splice[SN].site_score = splice_socre(1,1,chr,CT[k]-5)+splice_socre(1,2,chr,AC[j]-5);
            if(splice[SN].end-opt->change_length<=splice[SN].start) continue;
            if(SN>=200)continue;

            length_r = (splice[SN].start-start)+(end+1-splice[SN].end);
            if(length_r>=MAX_READ_LENGTH) continue;
            for(i = start;i<splice[SN].start;i++)
                ref[i-start] = chr[i];
            istart = i-start;
            for(i = splice[SN].end;i<end+1;i++)
                ref[istart+i-splice[SN].end] = chr[i];
            ref[length_r] = '\0';

            if(length_r==0)
            {
                splice[SN].cigar_num = 1;
                splice[SN].cigar = (struct cigar_t *)malloc(1*sizeof(struct cigar_t));
                splice[SN].cigar[0].c ='I';
                splice[SN].cigar[0].l =length;
                splice[SN].score = -4-(splice[SN].cigar[0].l-1)*opt->gap;
                SN++;
                if(SN>=200){(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free_SP(SN,splice);free(t_cigar);return 0;}
            }
            else
            {
                t_num = 0;
                splice[SN].score = NW(ref,length_r,text,length,t_cigar,&t_num,DP);
                splice[SN].cigar = (struct cigar_t *)malloc(t_num*sizeof(struct cigar_t));
                memcpy(splice[SN].cigar,t_cigar,t_num*sizeof(struct cigar_t));
                splice[SN].cigar_num = t_num;
                SN++;
                if(SN>=200){(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free_SP(SN,splice);free(t_cigar);return 0;}
            }
        }
    }
    }

    int best = 0;
    int best_score = splice[0].score;
    for(i = 1;i<SN;i++)
    {
        if(splice[i].score>best_score)
        {
            best_score = splice[i].score;
            best = i;
        }
        //else if(splice[i].score==best_score)
        //{
            //if(splice[i].site_score>splice[best].site_score)
            //{
                //best_score = splice[i].score;
                //best = i;
            //}
        //}
    }
    if(SN!=0)
    {
        (*score) = splice[best].score;
        uint64_t ref_pos = start;
        int read_pos = 0;

        for(i = 0;i<splice[best].cigar_num;i++)
        {
            if((splice[best].cigar[i].c=='M')||(splice[best].cigar[i].c=='X')||(splice[best].cigar[i].c=='D'))
            {
                if(ref_pos+splice[best].cigar[i].l<=splice[best].start)
                {
                    cigar[(*cigar_num)].c = splice[best].cigar[i].c;
                    cigar[(*cigar_num)].l = splice[best].cigar[i].l;
                    (*cigar_num)++;
                    if((*cigar_num)>=max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free_SP(SN,splice);free(t_cigar);return 0;}

                    ref_pos+=splice[best].cigar[i].l;
                    read_pos+=splice[best].cigar[i].l;
                }
                else
                {
                    if(splice[best].start-ref_pos>0)
                    {
                        cigar[(*cigar_num)].c = splice[best].cigar[i].c;
                        cigar[(*cigar_num)].l = splice[best].start-ref_pos;
                        (*cigar_num)++;
                        if((*cigar_num)>=max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free_SP(SN,splice);free(t_cigar);return 0;}
                        splice[best].cigar[i].l-= splice[best].start-ref_pos;
                    }
                    break;
                }
            }
            else if(splice[best].cigar[i].c=='I')
            {
                cigar[(*cigar_num)].c = splice[best].cigar[i].c;
                cigar[(*cigar_num)].l = splice[best].cigar[i].l;
                (*cigar_num)++;
                if((*cigar_num)>=max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free_SP(SN,splice);free(t_cigar);return 0;}
            }
            else {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free_SP(SN,splice);free(t_cigar);return 0;}
        }
        cigar[(*cigar_num)].c = 'N';
        cigar[(*cigar_num)].l = splice[best].length;
        (*cigar_num)++;
        if((*cigar_num)>=max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free_SP(SN,splice);free(t_cigar);return 0;}
        for(j = i;j<splice[best].cigar_num;j++)
        {
            cigar[(*cigar_num)].c = splice[best].cigar[j].c;
            cigar[(*cigar_num)].l = splice[best].cigar[j].l;
            (*cigar_num)++;
            if((*cigar_num)>=max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free_SP(SN,splice);free(t_cigar);return 0;}
        }
    }
    else {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free_SP(SN,splice);free(t_cigar);return 0;}

    free(ref);
    free(text);
    free(t_cigar);
    free_SP(SN,splice);

    return 0;
}
int Known_splice_pos1(int chr_order,uint64_t start,uint64_t end,char *seq,int length,struct cigar_t *cigar,int *cigar_num,int *score,int max,struct Splice_DP_t *DP)//MODE:1 GT-AG 2 CT-AC
{
    unsigned int i = 0,j;
    unsigned int Sstart,Send;

    int TH = 5;

    Sstart = start+length+TH;
    Send = end-length-TH;

    char *chr = opt->chr->list[chr_order].seq;
    unsigned int  chr_start= opt->chr->list[chr_order].start_site;

    char *text = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
    char *ref = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
    for(j = 0;j<length;j++)
        text[j] = seq[j];
    text[length] = '\0';

    struct cigar_t *t_cigar = (struct cigar_t *)calloc(MAX_CIGAR_BUF,sizeof(struct cigar_t));
    int t_num = 0;

    int istart;
    int length_r;

    struct splice_t splice[200];
    int SN = 0;
    int k = 0;

    int order = find_SP_Order(opt->SP,start+chr_start);
    if(order==-1) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free_SP(SN,splice);return 0;}

    for(k = order;k<opt->SP_num;k++)
    {
        if(opt->SP[k].start>Sstart+chr_start) break;
        if(opt->SP[k].end<Send+chr_start) continue;
        if(opt->SP[k].end>end+chr_start) continue;

        {
            splice[SN].start = opt->SP[k].start-chr_start;
            splice[SN].end = opt->SP[k].end-chr_start;
            splice[SN].length = splice[SN].end-splice[SN].start;
            //splice[SN].site_score = splice_socre(1,1,chr,CT[k]-5)+splice_socre(1,2,chr,AC[j]-5);
            if(SN>=200)break;

            length_r = (splice[SN].start-start)+(end+1-splice[SN].end);
            if(length_r>=MAX_READ_LENGTH) continue;
            ref[length_r] = '\0';
            for(i = start;i<splice[SN].start;i++)
                ref[i-start] = chr[i];
            istart = i-start;
            for(i = splice[SN].end;i<end+1;i++)
                ref[istart+i-splice[SN].end] = chr[i];

            if(length_r==0)
            {
                splice[SN].cigar_num = 1;
                splice[SN].cigar = (struct cigar_t *)malloc(1*sizeof(struct cigar_t));
                splice[SN].cigar[0].c ='I';
                splice[SN].cigar[0].l =length;
                splice[SN].score = -4-(splice[SN].cigar[0].l-1)*opt->gap;
                SN++;
                if(SN>=200){(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free_SP(SN,splice);return 0;}
            }
            else
            {
                t_num = 0;
                splice[SN].score = NW(ref,length_r,text,length,t_cigar,&(t_num),DP);
                splice[SN].score+=4+opt->SP[k].num;
                splice[SN].cigar = (struct cigar_t *)malloc(t_num*sizeof(struct cigar_t));
                memcpy(splice[SN].cigar,t_cigar,t_num*sizeof(struct cigar_t));
                splice[SN].cigar_num = t_num;
                SN++;
                if(SN>=200){(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free_SP(SN,splice);return 0;}
            }
        }
    }

    int best = 0;
    int best_score = splice[0].score;
    for(i = 1;i<SN;i++)
    {
        if(splice[i].score>best_score)
        {
            best_score = splice[i].score;
            best = i;
        }
        //else if(splice[i].score==best_score)
        //{
            //if(splice[i].site_score>splice[best].site_score)
            //{
                //best_score = splice[i].score;
                //best = i;
            //}
        //}
    }
    if(SN!=0)
    {
        (*score) = splice[best].score;
        uint64_t ref_pos = start;
        int read_pos = 0;

        for(i = 0;i<splice[best].cigar_num;i++)
        {
            if((splice[best].cigar[i].c=='M')||(splice[best].cigar[i].c=='X')||(splice[best].cigar[i].c=='D'))
            {
                if(ref_pos+splice[best].cigar[i].l<=splice[best].start)
                {
                    cigar[(*cigar_num)].c = splice[best].cigar[i].c;
                    cigar[(*cigar_num)].l = splice[best].cigar[i].l;
                    (*cigar_num)++;
                    if((*cigar_num)>max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free_SP(SN,splice);return 0;}

                    ref_pos+=splice[best].cigar[i].l;
                    read_pos+=splice[best].cigar[i].l;
                }
                else
                {
                    if(splice[best].start-ref_pos>0)
                    {
                        cigar[(*cigar_num)].c = splice[best].cigar[i].c;
                        cigar[(*cigar_num)].l = splice[best].start-ref_pos;
                        (*cigar_num)++;
                        if((*cigar_num)>max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free_SP(SN,splice);return 0;}
                        splice[best].cigar[i].l-= splice[best].start-ref_pos;
                    }
                    break;
                }
            }
            else if(splice[best].cigar[i].c=='I')
            {
                cigar[(*cigar_num)].c = splice[best].cigar[i].c;
                cigar[(*cigar_num)].l = splice[best].cigar[i].l;
                (*cigar_num)++;
                if((*cigar_num)>max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free_SP(SN,splice);return 0;}
            }
            else {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free_SP(SN,splice);return 0;}
        }
        cigar[(*cigar_num)].c = 'N';
        cigar[(*cigar_num)].l = splice[best].length;
        (*cigar_num)++;
        if((*cigar_num)>max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free_SP(SN,splice);return 0;}
        for(j = i;j<splice[best].cigar_num;j++)
        {
            cigar[(*cigar_num)].c = splice[best].cigar[j].c;
            cigar[(*cigar_num)].l = splice[best].cigar[j].l;
            (*cigar_num)++;
            if((*cigar_num)>max) {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free_SP(SN,splice);return 0;}
        }
    }
    else {(*cigar_num) = 0;(*score) = -1000;free(ref);free(text);free(t_cigar);free_SP(SN,splice);return 0;}

    free(ref);
    free(text);
    free(t_cigar);
    free_SP(SN,splice);
    return 0;
}

void extend_seed_tail_forward(int chr_order,char *seq,int length,int mode,unsigned int start,unsigned int end,struct seed_t *seed,int *seed_num,unsigned int pos,int Sstart,int init_N)
{
    struct seed_t seed_r[10];
    int flag = 0;

    //if(Sstart>=5)
    if(Sstart>0)
        flag = tail_seed(4,chr_order,0,pos+opt->chr->list[chr_order].start_site,seq,Sstart,0,seed_r);
    if(flag)
    {
        int i;
        for(i = 0;i<flag;i++)
        {
        seed_r[i].pos -=opt->chr->list[chr_order].start_site;
        seed_r[i].abs = seed_r[i].pos-seed_r[i].start;
        if((mode!=1)&&(seed_r[i].pos<start))
        {
            seed_r[i].start += start-seed_r[i].pos+1;
            seed_r[i].length -= start-seed_r[i].pos+1;
            seed_r[i].pos += start-seed_r[i].pos+1;
        }
        if(seed_r[i].length>0)
        {
            insert_seed_start(seed,seed_num,init_N,SEED_BUF_LENGTH,seed_r[i]);

        if(i==0)
            extend_seed_tail_forward(chr_order,seq,seed_r[i].start,mode,start,seed_r[i].pos,seed,seed_num,seed_r[i].pos,seed_r[i].start,init_N);
        }
        }
    }
}
void extend_seed_tail_backward(int chr_order,char *seq,int length,int mode,unsigned int start,unsigned int end,struct seed_t *seed,int *seed_num,unsigned int pos,int Sstart,int init_N)
{
    struct seed_t seed_r[10];
    int flag = 0;

    //if(length-read_pos>=5)
    if(length-Sstart>0)
        flag = tail_seed(4,chr_order,0,pos-1+opt->chr->list[chr_order].start_site,&seq[Sstart],length-Sstart,1,seed_r);
    if(flag)
    {
        int i;
        for(i = 0;i<flag;i++)
        {
            seed_r[i].pos -= opt->chr->list[chr_order].start_site;
        seed_r[i].start += Sstart;
        seed_r[i].abs = seed_r[i].pos-seed_r[i].start;
        if((mode!=0)&&(seed_r[i].pos+seed_r[i].length>end))
            seed_r[i].length-=seed_r[i].pos+seed_r[i].length-end;
        if(seed_r[i].length>0)
        {
            insert_seed_start(seed,seed_num,init_N,SEED_BUF_LENGTH,seed_r[i]);

        if(i==0)
            extend_seed_tail_backward(chr_order,seq,length,mode,seed_r[i].pos+seed_r[i].length-1,end,seed,seed_num,seed_r[i].pos+seed_r[i].length,seed_r[i].start+seed_r[i].length,init_N);
        }

        }
    }
}
void extend_seed_tail(char *chr,int chr_order,char *seq,int mode,unsigned int start,unsigned int end,struct seed_t *seed,int *seed_num,int init_N)
{
    int i;
    int init_num = (*seed_num);
    int length = strlen(seq);

    for(i = 0;i<init_num;i++)
    {
        if((seed[i].start == 0)&&(seed[i].pos == 0)&&(seed[i].length == 0)) continue;
        extend_seed_tail_forward(chr_order,seq,length,mode,start,seed[i].pos,seed,seed_num,seed[i].pos,seed[i].start,init_N);
        extend_seed_tail_backward(chr_order,seq,length,mode,seed[i].pos+seed[i].length,end,seed,seed_num,seed[i].pos+seed[i].length,seed[i].start+seed[i].length,init_N);
    }
    //if(mode!=0) extend_seed_tail_forward(chr_order,seq,length,mode,start,end,seed,seed_num,end,length);
    //if(mode!=1) extend_seed_tail_backward(chr_order,seq,length,mode,start,end,seed,seed_num,start,0);
}
void cigar2site(unsigned int pos,unsigned int *site,struct cigar_t *cigar,int cigar_num)
{
    int read_site = 0;
    unsigned int ref_site = pos;
    int i,j;

    for(i = 0;i<cigar_num;i++)
    {
        switch(cigar[i].c)
        {
        case 'M':
            for(j = 0;j<cigar[i].l;j++)
                site[read_site+j] = ref_site+j;

            ref_site+=cigar[i].l;
            read_site+=cigar[i].l;
            break;
        case 'X':
            for(j = 0;j<cigar[i].l;j++)
                site[read_site+j] = ref_site+j;

            ref_site+=cigar[i].l;
            read_site+=cigar[i].l;
            break;
        case 'I':
            for(j = 0;j<cigar[i].l;j++)
                site[read_site+j] = ref_site-1;

            read_site+=cigar[i].l;
            break;
        case 'D':
            ref_site+=cigar[i].l;
            break;
        case 'N':
            ref_site+=cigar[i].l;
            break;
        case 'U':
            ref_site+=cigar[i].l;
            break;
        case 'S':
            for(j = 0;j<cigar[i].l;j++)
                site[read_site+j] = 0;
            read_site+=cigar[i].l;
            break;
        }
    }
}
void cigar2site2(unsigned int pos,unsigned int *site,int length,struct cigar_t *cigar,int cigar_num)
{
    int read_site = length-1;
    unsigned int ref_site = pos;
    int i,j;

    for(i = cigar_num-1;i>=0;i--)
    {
        switch(cigar[i].c)
        {
        case 'M':
        case 'X':
            for(j = 0;j<cigar[i].l;j++)
                site[read_site-j] = ref_site-j;

            ref_site-=cigar[i].l;
            read_site-=cigar[i].l;
            break;
        case 'I':
            for(j = 0;j<cigar[i].l;j++)
                site[read_site-j] = ref_site;

            read_site-=cigar[i].l;
            break;
        case 'D':
            ref_site-=cigar[i].l;
            break;
        case 'N':
            ref_site-=cigar[i].l;
            break;
        case 'U':
            ref_site-=cigar[i].l;
            break;
        case 'S':
            for(j = 0;j<cigar[i].l;j++)
                site[read_site-j] = 0;
            read_site-=cigar[i].l;
            break;
        }
    }
}
#define AREA_SEED_CAND_NUM 30
#define AREA 500000
int cigar_score(struct cigar_t *cigar,int cigar_num)
{
    int i;
    int x = cigar_num;
    int score = 0;
    for(i = 0;i<x;i++)
    {
        if(cigar[i].c=='M') score+=cigar[i].l*opt->match;
        else if(cigar[i].c=='X') score-=cigar[i].l*opt->miss;
        else if(cigar[i].c=='I') score-=4+(cigar[i].l-1)*opt->insert;
        else if(cigar[i].c=='D')
        {
            if(cigar[i].l>opt->change_length){cigar[i].c='N';score-=opt->splice;}
            else score-=4+(cigar[i].l-1)*opt->gap;
        }
        else if(cigar[i].c=='N') score-=opt->splice;
    }
	return score;
}

void check_cigar(char *read,char *cigarS,int length)
{
    struct cigar_t *cigar = (struct cigar_t *)calloc(MAX_CIGAR_BUF,sizeof(struct cigar_t));
    int cigar_num = 0;

    char temp[20];
    int l = strlen(cigarS);
    int temp_length = 0;
    int i = 0,j = 0;

    cigar[cigar_num].l = 0;
    j = 0;
    for(i = 0;i<l;i++)
    {
        if((cigarS[i]>='0')&&(cigarS[i]<='9'))
        {
            temp[j] = cigarS[i];
            j++;
            if(j>=20)
                printf("WRONG CIGAR LENGTH\n");
        }
        else
        {
            temp[j] = '\0';
            cigar[cigar_num].c = cigarS[i];
            cigar[cigar_num].l = atoi(temp);
            if(cigar[cigar_num].l<=0)
                printf("WRONG CIGAR LENGTH\n");
            cigar_num++;
            cigar[cigar_num].l = 0;
            j = 0;
        }
    }
    temp_length=0;
    for(i = 0;i<cigar_num;i++)
    {
        if(cigar[i].c =='M') temp_length+=cigar[i].l;
        //else if(cigar[i].c =='X') temp_length+=cigar[i].l;
        else if(cigar[i].c =='D') temp_length+=0;
        else if(cigar[i].c =='I')
        {
            temp_length+=cigar[i].l;
            //if(cigar[i].l>20)
                //printf("SUPER LONG INSERT\n");
        }
        else if(cigar[i].c =='N') temp_length+=0;
        else if(cigar[i].c =='S') temp_length+=cigar[i].l;
        else
            printf("WRONG CIGAR\n");
    }
    if(temp_length!= length)
        printf("WRONG READ LENGTH\n");
}
void check_cigar2(struct cigar_t *cigar,int cigar_num,int length)
{
    int temp_length = 0;
    int i = 0;

    temp_length=0;
    for(i = 0;i<cigar_num;i++)
    {
        if(cigar[i].l<=0)
            printf("WRONG CIGAR LENGTH\n");
        if(cigar[i].c =='M') temp_length+=cigar[i].l;
        else if(cigar[i].c =='X') temp_length+=cigar[i].l;
        else if(cigar[i].c =='D') temp_length+=0;
        else if(cigar[i].c =='I')
        {
            temp_length+=cigar[i].l;
            //if(cigar[i].l>20)
                //printf("SUPER LONG INSERT\n");
        }
        else if(cigar[i].c =='N') temp_length+=0;
        else if(cigar[i].c =='S') temp_length+=cigar[i].l;
        else
            printf("WRONG CIGAR\n");
    }
    if(temp_length!= length)
        printf("WRONG READ LENGTH\n");
}
void Update_SP_s(struct splice_list *SP,int *SP_num,unsigned int MAX_N,struct splice_list *sp,int flag)
{
    if(((*SP_num)>=MAX_N)||(flag==1))
    {
        pthread_mutex_lock(&UpdateLock);
        qsort(SP,(*SP_num),sizeof(struct splice_list),SP_cmp);
        int i,j = 0;
        int SP_total = opt->SP_num;
        for(i = 0;i<(*SP_num);i++)
        {
            flag = 0;
            while(j<SP_total)
            {
                if(opt->SP[j].start>SP[i].start) break;
                if((opt->SP[j].start==SP[i].start)&&(opt->SP[j].end==SP[i].end))
                    {opt->SP[j].num++;flag = 1;break;}
                j++;
            }
            if(flag==0) Update_SP(SP+i);
            while((opt->SP[j].start>=SP[i].start)&&(j>0))j--;
        }
        qsort(opt->SP,opt->SP_num,sizeof(struct splice_list),SP_cmp);
        (*SP_num) = 0;
        pthread_mutex_unlock(&UpdateLock);
    }
    memcpy(SP+(*SP_num),sp,1*sizeof(struct splice_list));
    (*SP_num) ++;
}
void check_exon(struct read_t *read0,char *name,char *read,int length,struct cigar_t *cigar,int *cigar_num,uint64_t *start_pos,int chr_order,struct Splice_DP_t *DP,char *cigarS,struct splice_list *SP,int *SP_num)
{
    if(strcmp(name,"SimG2_S2294_4")==0)
        printf("1\n");
    int break_flag = 1;
    int EXON_LENGTH = 50;

    int i,j;
    unsigned int ref_pos = *start_pos;
    int read_pos = 0;

    struct exon_seed_t *seed = (struct exon_seed_t *)malloc(EXON_SEED*sizeof(struct exon_seed_t));
    int seed_num = 0;
    int n = 0;
    int score = 0;

    struct cigar_t *t_cigar = (struct cigar_t *)malloc(MAX_CIGAR_BUF*sizeof(struct cigar_t));
    int t_num = 0;

    int seed_order[1000];
    int seed_order_n = 0;

    int max_s,max_o;

    seed[0].ref_start = *start_pos;
    seed[0].read_start = 0;
    seed[0].flag = 0;
    seed[0].front = -1;
    for(i = 0;i<(*cigar_num);i++)
    {
        if (cigar[i].c=='M') {ref_pos+=cigar[i].l;read_pos+=cigar[i].l;n++;score+=opt->match*cigar[i].l;}
        else if(cigar[i].c=='X'){ref_pos+=cigar[i].l;read_pos+=cigar[i].l;n++;score-=opt->miss*cigar[i].l;}
        else if(cigar[i].c=='I'){read_pos+=cigar[i].l;n++;score-=4+(cigar[i].l-1)*opt->gap;}
        else if(cigar[i].c=='S'){read_pos+=cigar[i].l;n++;}
        else if(cigar[i].c=='D'){ref_pos+=cigar[i].l;n++;score-=4+(cigar[i].l-1)*opt->gap;}
        else if(cigar[i].c=='N')
        {
            seed[seed_num].ref_end = ref_pos-1;
            seed[seed_num].read_end = read_pos-1;
            seed[seed_num].score = score;
            seed[seed_num].cigar_num = n;
            if(seed[seed_num].ref_end-seed[seed_num].ref_start+1>=EXON_LENGTH)seed[seed_num].flag = 1;
            seed[seed_num].cigar = (struct cigar_t *)malloc(n*sizeof(struct cigar_t));
            memcpy(seed[seed_num].cigar,(cigar+i-n),n*sizeof(struct cigar_t));
            break_flag = break_flag&seed[seed_num].flag;
            seed_num++;

            ref_pos+=cigar[i].l;
            score = 0;n = 0;
            seed[seed_num].ref_start = ref_pos;
            seed[seed_num].read_start = read_pos;
            seed[seed_num].flag = 0;
            seed[seed_num].front = seed_num-1;
        }
    }
    {
    seed[seed_num].ref_end = ref_pos-1;
    seed[seed_num].read_end = read_pos-1;
    seed[seed_num].score = score;
    seed[seed_num].cigar_num = n;
    if(seed[seed_num].ref_end-seed[seed_num].ref_start+1>=EXON_LENGTH)seed[seed_num].flag = 1;
    seed[seed_num].cigar = (struct cigar_t *)malloc(n*sizeof(struct cigar_t));
    memcpy(seed[seed_num].cigar,(cigar+i-n),n*sizeof(struct cigar_t));
    break_flag = break_flag&seed[seed_num].flag;
    seed_num++;
    }

    if(read_pos!=length)
        printf("1\n");

    if(!break_flag)
    {
    struct cigar_t * temp_cigar;

    int init_num = seed_num;
    int front = 0;
    int back = init_num-1;

    i = 1;
    while(i<init_num)
    {
        if(seed[i].flag!=1)
        {
            if((seed[i].read_end+1-seed[i].read_start)==0)
            {
                if(i+1>=init_num) back = front;
                else seed[i+1].front = front;
                i++;
                continue;
            }
            seed[seed_num].ref_start = seed[front].ref_end+1;
            seed[seed_num].read_start = seed[i].read_start;
            seed[seed_num].ref_end = 0;
            seed[seed_num].read_end = seed[i].read_end;
            seed[seed_num].score = 0;

            t_num = 0;
            if(i+1<init_num)
            seed[seed_num].score = Tail_DP(opt->chr->list[chr_order].seq,seed[seed_num].ref_start,seed[i+1].ref_start,1,
                                           &(read[seed[seed_num].read_start]),(seed[seed_num].read_end-seed[seed_num].read_start+1),
                                               t_cigar,&t_num,DP,&(seed[seed_num].ref_end),0);
            else
                seed[seed_num].score = Tail_DP(opt->chr->list[chr_order].seq,seed[seed_num].ref_start,opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length,1,
                                               &(read[seed[seed_num].read_start]),(seed[seed_num].read_end-seed[seed_num].read_start+1),
                                               t_cigar,&t_num,DP,&(seed[seed_num].ref_end),1);

            if((seed[seed_num].score>=seed[i].score-opt->splice)||((i == init_num-1)&&(seed[i].ref_end-seed[i].ref_start<opt->cut_tail)))
            {
                if(t_cigar[0].c==seed[front].cigar[seed[front].cigar_num-1].c)
                {
                    seed[front].cigar[seed[front].cigar_num-1].l+= t_cigar[0].l;
                    temp_cigar = (struct cigar_t *)malloc((t_num-1+seed[front].cigar_num)*sizeof(struct cigar_t));
                    memcpy(temp_cigar,seed[front].cigar,seed[front].cigar_num*sizeof(struct cigar_t));
                    memcpy(temp_cigar+seed[front].cigar_num,t_cigar+1,(t_num-1)*sizeof(struct cigar_t));
                    seed[front].cigar_num+=t_num-1;
                }
                else
                {
                    temp_cigar = (struct cigar_t *)malloc((t_num+seed[front].cigar_num)*sizeof(struct cigar_t));
                    memcpy(temp_cigar,seed[front].cigar,seed[front].cigar_num*sizeof(struct cigar_t));
                    memcpy(temp_cigar+seed[front].cigar_num,t_cigar,t_num*sizeof(struct cigar_t));
                    seed[front].cigar_num+=t_num;
                }
                free(seed[front].cigar);
                seed[front].cigar = temp_cigar;

                seed[front].ref_end = seed[seed_num].ref_end;
                seed[front].read_end = seed[seed_num].read_end;
                seed[front].score+=seed[seed_num].score;

                if(i+1>=init_num) back = front;
                else seed[i+1].front = front;
            }
            else front = i;
        }
        else front = i;
        i++;
    }

    i = max_s= max_o = back;
    while(i!=-1)
    {
        back = i;
        i = seed[back].front;
        if((i!=-1)&&(seed[i].flag!=1))
        {
            if((seed[i].read_end+1-seed[i].read_start)==0)
            {
                seed[back].front = seed[i].front;
                i = back;
                continue;
            }
            seed[seed_num].ref_start = 0;
            seed[seed_num].read_start = seed[i].read_start;
            seed[seed_num].ref_end = seed[back].ref_start-1;
            seed[seed_num].read_end = seed[i].read_end;
            seed[seed_num].score = 0;

            t_num = 0;
            if(seed[i].front!=-1)
            seed[seed_num].score = Tail_DP(opt->chr->list[chr_order].seq,seed[seed_num].ref_end,seed[seed[i].front].ref_end,0,
                                           &(read[seed[seed_num].read_start]),(seed[seed_num].read_end-seed[seed_num].read_start+1)
                                       ,t_cigar,&t_num,DP,&(seed[seed_num].ref_start),0);
            else
            seed[seed_num].score = Tail_DP(opt->chr->list[chr_order].seq,seed[seed_num].ref_end,0,0,
                                           &(read[seed[seed_num].read_start]),(seed[seed_num].read_end-seed[seed_num].read_start+1)
                                       ,t_cigar,&t_num,DP,&(seed[seed_num].ref_start),1);

            if((seed[seed_num].score>=seed[i].score-opt->splice)||((i==0)&&(seed[i].ref_end-seed[i].ref_start<opt->cut_tail)))
            {
                if(t_cigar[t_num-1].c==seed[back].cigar[0].c)
                {
                    t_cigar[t_num-1].l+=seed[back].cigar[0].l;
                    temp_cigar = (struct cigar_t *)malloc((t_num-1+seed[back].cigar_num)*sizeof(struct cigar_t));
                    memcpy(temp_cigar,t_cigar,t_num*sizeof(struct cigar_t));
                    memcpy(temp_cigar+t_num,seed[back].cigar+1,(seed[back].cigar_num-1)*sizeof(struct cigar_t));
                    seed[back].cigar_num+=t_num-1;
                }
                else
                {
                    temp_cigar = (struct cigar_t *)malloc((t_num+seed[back].cigar_num)*sizeof(struct cigar_t));
                    memcpy(temp_cigar,t_cigar,t_num*sizeof(struct cigar_t));
                    memcpy(temp_cigar+t_num,seed[back].cigar,(seed[back].cigar_num)*sizeof(struct cigar_t));
                    seed[back].cigar_num+=t_num;
                }
                free(seed[back].cigar);
                seed[back].cigar = temp_cigar;

                seed[back].ref_start = seed[seed_num].ref_start;
                seed[back].read_start = seed[seed_num].read_start;
                seed[back].score+=seed[seed_num].score;
                seed[back].front = seed[i].front;
                i = back;
            }
        }
    }

    i = back= max_s = max_o;
    while(i!=-1)
    {
        back = i;
        i = seed[back].front;
        if((i!=-1)&&(seed[i].front!=-1)&&(seed[i].flag!=1)
           &&(seed[back].ref_start-seed[i].ref_end>opt->change_length+1)
           &&(seed[i].ref_start-seed[seed[i].front].ref_end>opt->change_length+1))
        {
            seed[seed_num].ref_start = seed[seed[i].front].ref_end+1;
            seed[seed_num].read_start = seed[seed[i].front].read_end+1;
            seed[seed_num].ref_end = seed[back].ref_start-1;
            seed[seed_num].read_end = seed[back].read_start-1;
            seed[seed_num].score = 0;

            t_num = 0;
            seed[seed_num].score = Splice_DP(opt->chr->list[chr_order].seq,seed[seed_num].ref_start,seed[seed_num].ref_end,
                            &(read[seed[seed_num].read_start]),(seed[seed_num].read_end-seed[seed_num].read_start+1),t_cigar,&t_num,DP);
            if(seed[seed_num].score>seed[i].score-2*opt->splice)
            {
                read_pos = seed[seed_num].read_start;
                ref_pos = seed[seed_num].ref_start;
                for(j = 0;j<t_num;j++)
                {
                    if((t_cigar[j].c=='M')||(t_cigar[j].c=='X')) {read_pos+=t_cigar[j].l;ref_pos+=t_cigar[j].l;}
                    else if((t_cigar[j].c=='I')||(t_cigar[j].c=='S')) {read_pos+=t_cigar[j].l;}
                    else if(t_cigar[j].c=='D') {ref_pos+=t_cigar[j].l;}
                    if(t_cigar[j].c=='N') break;
                }
                if(seed[seed[i].front].cigar[seed[seed[i].front].cigar_num-1].c==t_cigar[0].c)
                {
                    temp_cigar = (struct cigar_t *)malloc((j-1+seed[seed[i].front].cigar_num)*sizeof(struct cigar_t));
                    memcpy(temp_cigar,seed[seed[i].front].cigar,seed[seed[i].front].cigar_num*sizeof(struct cigar_t));
                    temp_cigar[seed[seed[i].front].cigar_num-1].l+=t_cigar[0].l;
                    memcpy(temp_cigar+seed[seed[i].front].cigar_num,t_cigar+1,(j-1)*sizeof(struct cigar_t));
                    seed[seed[i].front].cigar_num += j-1;
                    free(seed[seed[i].front].cigar);
                    seed[seed[i].front].cigar = temp_cigar;
                }
                else
                {
                    temp_cigar = (struct cigar_t *)malloc((j+seed[seed[i].front].cigar_num)*sizeof(struct cigar_t));
                    memcpy(temp_cigar,seed[seed[i].front].cigar,seed[seed[i].front].cigar_num*sizeof(struct cigar_t));
                    memcpy(temp_cigar+seed[seed[i].front].cigar_num,t_cigar,j*sizeof(struct cigar_t));
                    seed[seed[i].front].cigar_num += j;
                    free(seed[seed[i].front].cigar);
                    seed[seed[i].front].cigar = temp_cigar;
                }
                seed[seed[i].front].ref_end = ref_pos-1;
                seed[seed[i].front].read_end = read_pos-1;

                ref_pos+=t_cigar[j].l;

                seed[back].ref_start = ref_pos;
                seed[back].read_start = read_pos;

                j++;
                if(t_cigar[t_num-1].c==seed[back].cigar[0].c)
                {
                    temp_cigar = (struct cigar_t *)malloc((t_num-j-1+seed[back].cigar_num)*sizeof(struct cigar_t));
                    memcpy(temp_cigar,t_cigar+j,(t_num-j)*sizeof(struct cigar_t));
                    temp_cigar[t_num-j-1].l+=seed[back].cigar[0].l;
                    memcpy(temp_cigar+t_num-j,seed[back].cigar+1,(seed[back].cigar_num-1)*sizeof(struct cigar_t));
                    seed[back].cigar_num += t_num-j-1;
                    free(seed[back].cigar);
                    seed[back].cigar = temp_cigar;
                }
                else
                {
                    temp_cigar = (struct cigar_t *)malloc((t_num-j+seed[back].cigar_num)*sizeof(struct cigar_t));
                    memcpy(temp_cigar,t_cigar+j,(t_num-j)*sizeof(struct cigar_t));
                    memcpy(temp_cigar+t_num-j,seed[back].cigar,(seed[back].cigar_num)*sizeof(struct cigar_t));
                    seed[back].cigar_num += t_num-j;
                    free(seed[back].cigar);
                    seed[back].cigar = temp_cigar;
                }

                seed[back].front = seed[i].front;
                i = back;
            }
        }
    }

    max_s = max_o;
    }
    else max_o = seed_num-1;

    seed_order_n = 0;
    while(max_o!=-1)
    {
        seed_order[seed_order_n] = max_o;
        seed_order_n++;
        max_o = seed[max_o].front;
    }
    for(i = 0;i<seed_order_n/2;i++)
    {
        max_o = seed_order[i];
        seed_order[i]=seed_order[seed_order_n-1-i];
        seed_order[seed_order_n-1-i] = max_o;
    }
    if(seed_order[0]!=0)
        *start_pos = seed[seed_order[0]].ref_start;

    //write cigar
    char temp[20];
    cigarS[0] = '\0';

    unsigned int start = 0;
    struct splice_list sp;
    unsigned int chr_start = opt->chr->list[chr_order].start_site;

    for(i = 0;i<seed_order_n;i++)
    {
        start = seed[seed_order[i]].ref_start;
        for(j = 0;j<seed[seed_order[i]].cigar_num;j++)
        {
            if((seed[seed_order[i]].cigar[j].l<=8)&&(seed[seed_order[i]].cigar[j].c=='N')) seed[seed_order[i]].cigar[j].c = 'D';
            if((seed[seed_order[i]].cigar[j].l>8)&&(seed[seed_order[i]].cigar[j].c=='D')) seed[seed_order[i]].cigar[j].c = 'N';
            sprintf(temp,"%d%c",seed[seed_order[i]].cigar[j].l,seed[seed_order[i]].cigar[j].c);
            strcat(cigarS,temp);

            if((seed[seed_order[i]].cigar[j].c=='M')||(seed[seed_order[i]].cigar[j].c=='X')) start+=seed[seed_order[i]].cigar[j].l;
            if(seed[seed_order[i]].cigar[j].c == 'N')
            {
                sp.start = start+chr_start;
                sp.end = start+seed[seed_order[i]].cigar[j].l+chr_start;
                sp.num = 1;
                Update_SP_s(SP,SP_num,opt->SP_block_MAX,&sp,0);

                start+=seed[seed_order[i]].cigar[j].l;
            }
            else if(seed[seed_order[i]].cigar[j].c == 'D') start+=seed[seed_order[i]].cigar[j].l;
        }
        if(i+1<seed_order_n)
        {
            sprintf(temp,"%dN",seed[seed_order[i+1]].ref_start-seed[seed_order[i]].ref_end-1);
            strcat(cigarS,temp);
        }
    }

    for(i = 0;i<seed_num;i++)
    {
        if(seed[i].cigar!=NULL)free(seed[i].cigar);
    }
    free(seed);
    free(t_cigar);
}

int find8(char *chr,char *seq,int mode,unsigned int start,unsigned int end,
               unsigned int **hash,int *hash_num,struct seed_t *seed,unsigned int *site,int align_flag,int chr_order,struct Splice_DP_t *DP)
{
    int i,j,k;
    int length = strlen(seq);

    int length_t;

    char f_seq[KMER_FRONT+1];
    uint32_t f_base = 0;
    uint32_t fmax = pow(4,KMER_FRONT);
    int hash_site = 0;

    int seed_num = 0;

    int temp_start = 0;
    int temp_end = 0;
    int temp_mode = 0;
    unsigned int start_site = 0;
    unsigned int end_site = 0;
    int flag = 0;

    int step = 0;

    unsigned int chr_start = opt->chr->list[chr_order].start_site;
    unsigned int START = start+chr_start;
    unsigned int END = end+chr_start;
    unsigned int start_pos = 0;

    if(align_flag==1)
    {
        struct cigar_t *t_cigar = (struct cigar_t *)malloc((MAX_CIGAR_BUF+2)*sizeof(struct cigar_t));
        int t_num;

        if(mode==0)//((mode==0)&&(length<200))
        {
            Tail_DP(chr,start+1,end,1,seq,length,t_cigar,&t_num,DP,&start_pos,1);
            cigar2site(start+1,site,t_cigar,t_num);
        }
        else if(mode==1)//if((mode==1)&&(length<200))
        {
            Tail_DP(chr,end-1,start,0,seq,length,t_cigar,&t_num,DP,&start_pos,1);
            cigar2site2(end-1,site,length,t_cigar,t_num);
        }
        else if(mode==2)
        {
            if((end-start-1)<=2*max(length*(1.1),length+20))
                Area_DP(chr,start+1,end-1,seq,length,t_cigar,&t_num,DP);
            else
                Splice_DP(chr,start+1,end-1,seq,length,t_cigar,&t_num,DP);
            cigar2site(start+1,site,t_cigar,t_num);
        }
        free(t_cigar);
    }
    else
    {
    {
        step = 1;//KMER_FRONT/2;
        i = 0;
        while ((i+KMER_FRONT)<=length)
        {
            strncpy(f_seq,&seq[i],KMER_FRONT);
            f_seq[KMER_FRONT] = '\0';
            f_base = char2base_16(f_seq,KMER_FRONT);

            if(f_base>=fmax) {return -1;}

            hash_site = hash_c_find(hash,hash_num,f_base,START);
            if((hash_site+AREA_SEED_CAND_NUM<hash_num[f_base])&&(hash[f_base][hash_site+AREA_SEED_CAND_NUM]>=START)&&(hash[f_base][hash_site+AREA_SEED_CAND_NUM]<=END)) {i+=step;continue;}
            while ((hash_site < hash_num[f_base])&&(hash[f_base][hash_site]<=END))
            {
                if((hash[f_base][hash_site]>=START)&&(hash[f_base][hash_site]<=END))
                {
                    seed[seed_num].pos = hash[f_base][hash_site]-chr_start;
                    seed[seed_num].start = i;
                    seed[seed_num].length = KMER_FRONT;
                    seed[seed_num].abs = seed[seed_num].pos-seed[seed_num].start;
                    seed[seed_num].last = -1;
                    seed[seed_num].score = 0;
                    seed_num++;
                    //insert_seed(seed,&seed_num,SEED_BUF_LENGTH,seed_r);
                }
                hash_site++;
            }
            i+=step;
        }
    }


        int s = 0;
        int l = 0;
        //seed前后延伸
        for(i = 0;i<seed_num;i++)
        {
            k = 0;
            s = seed[i].start;
            for(j = 1;j<=s;j++)
            {
                if(nst_nt4_table[(int)chr[seed[i].pos-j]]!=nst_nt4_table[(int)seq[s-j]])
                    break;
                k = j;
            }
            seed[i].pos-=k;
            seed[i].start-=k;
            seed[i].length+=k;

            s = seed[i].start;
            l = seed[i].length;
            k = 0;
            for(j = 0;j<length-s-l;j++)
            {
                if(nst_nt4_table[(int)chr[seed[i].pos+j+l]]!=nst_nt4_table[(int)seq[s+l+j]])
                    break;
                k = j+1;
            }
            seed[i].length+=k;
        }
        qsort(seed,seed_num,sizeof(struct seed_t),seed_cmp);

        int seed_start = 0;
        for(i = 0;i<seed_num;i++)
        {
            while((seed[i].start >seed[seed_start].start)&&(seed_start<seed_num)) seed_start++;
            for(j = seed_start;j<i;j++)
            {
                if(seed[j].length == 0)continue;
                if((seed[i].start == seed[j].start)&&(seed[i].pos == seed[j].pos))
                {seed[i].start = 0;seed[i].pos = 0;seed[i].length = 0;seed[i].score = 0;seed[i].last = -1;break;}
            }
        }
        qsort(seed,seed_num,sizeof(struct seed_t),seed_cmp);
        int init_num = seed_num;
        extend_seed_tail(chr,chr_order,seq,mode,start,end,seed,&seed_num,init_num);

        int max_o = -1;
        int max_s = 0;

        if(mode!=1) qsort(seed,seed_num,sizeof(struct seed_t),seed_cmp);
        else qsort(seed,seed_num,sizeof(struct seed_t),seed_cmp_r);

        seed_start = 0;
        for(i = 0;i<seed_num;i++)
        {
            while((seed[i].start >seed[seed_start].start)&&(seed_start<seed_num)) seed_start++;
            for(j = seed_start;j<i;j++)
            {
                if(seed[j].length == 0)continue;
                if((seed[i].start == seed[j].start)&&(seed[i].pos == seed[j].pos))
                {seed[i].start = 0;seed[i].pos = 0;seed[i].length = 0;seed[i].score = 0;seed[i].last = -1;break;}
            }
        }

        for(i = 0; i<seed_num; i++)
        {
            if(seed[i].length == 0){seed[i].score = -1000;continue;}

            seed[i].score = seed[i].length*(opt->match);
            if(mode!=1)
            {
                if((seed[i].abs>=start+1)&&(seed[i].abs-start>opt->change_length)) seed[i].score-=opt->splice;
                //else if ((seed[i].abs<start+1)&&(start+1-seed[i].abs<15)) {seed[i].score -= 4+opt->gap*(start-1-seed[i].abs);}
                //else if ((seed[i].abs<start+1)&&(start+1-seed[i].abs>=15)) {seed[i].score = -1000;continue;}
                else if ((seed[i].abs>=start+1)&&(seed[i].abs-start+1<=opt->change_length)) seed[i].score-=4+(seed[i].abs-start-1-1)*opt->gap;
                else if ((seed[i].abs<start+1)) seed[i].score-=4+(start+1-seed[i].abs)*opt->match+(start+1-seed[i].abs-1)*opt->gap;
            }

            max_o = -1;
            max_s = 0;
            for(j = 0;j<i;j++)
            {
                if((seed[i].start > seed[j].start)&&((seed[i].start+seed[i].length)>(seed[j].start+seed[j].length)&&(seed[i].pos > seed[j].pos)))
                {
                    if(seed[i].start<=seed[j].start+seed[j].length)
                        max_s = seed[j].score + (seed[i].start+seed[i].length-seed[j].start-seed[j].length)*opt->match;
                    else
                        max_s = seed[j].score + seed[i].length*opt->match;

                    if(seed[i].abs<seed[j].abs)// max_s--;//max_s -=(seed[j].abs-seed[i].abs)+1;
                         max_s -=(seed[j].abs-seed[i].abs)*opt->match+4+(seed[j].abs-seed[i].abs-1)*opt->gap;
                    //{
                        //l = seed[j].abs-seed[i].abs;
                        //if(l>15) continue;
                        //else max_s -=(seed[j].abs-seed[i].abs)*opt->match+4+(seed[j].abs-seed[i].abs-1)*opt->gap;
                    //}

                    else if((seed[i].abs>seed[j].abs)&&(seed[i].abs-seed[j].abs<opt->change_length))
                        max_s -=4+(seed[i].abs-seed[j].abs-1)*opt->gap;
                    else if((seed[i].abs>seed[j].abs)&&(seed[i].abs-seed[j].abs>opt->change_length))
                        max_s -=opt->splice;

                    if(seed[i].score < max_s)
                    {
                        max_o = j;
                        seed[i].score = max_s;
                    }
                }
            }
            if(max_o != -1) seed[i].last = max_o;
            //if(max_o == -1)
            //{
                //seed[i].score = seed[i].length*(opt->match);
                //if((mode!=1)&&(seed[i].abs-start>opt->change_length)) seed[i].score-=opt->splice;
            //}
            //else
                //seed[i].last = max_o;

        }
        //回溯
        max_o = -1;
        max_s = -1000;
        int best_num = 0;
        for(i = seed_num-1; i>=0; i--)
        {
            if(mode!=0)
            {
                if((end-length>seed[i].abs)&&(end-length-seed[i].abs>opt->change_length)) seed[i].score-=opt->splice;
                //else if((end-length<seed[i].abs)&&(seed[i].abs-(end-length)>15)) continue;

                else if((end-length>seed[i].abs)&&(end-length-seed[i].abs<=opt->change_length)) seed[i].score-=4+(end-length-seed[i].abs-1)*opt->gap;
                else if((end-length<seed[i].abs)) seed[i].score-=4+(seed[i].abs-(end-length))*opt->match+(seed[i].abs-(end-length)-1)*opt->gap;
            }
            if(seed[i].score>=max_s)
            {
                max_o = i;
                max_s = seed[i].score;
                best_num = 1;
            }
        }
        if(best_num>1) max_o = -1;
        //if(max_s<0) max_o = -1;
        while(max_o!=-1)
        {
            for(i = 0;i<seed[max_o].length;i++)
                site[seed[max_o].start+i] = seed[max_o].pos + i;
            if(seed[max_o].last!=-1)
                max_o = seed[max_o].last;
            else max_o = -1;
        }
        for(j = 0;j<length;j++)
        {
            if(site[j]>=end) site[j] = end-1;
            else if((site[j]!=0)&&(site[j]<start)) site[j] = start;

            if((j!=0)&&(site[j]!=0)&&(site[j]<site[j-1]))
                site[j]=site[j-1];

            if ((flag==0)&&(site[j]==0))
            {
                temp_start = j;
                flag = 1;
            }
            if ((flag==1)&&((site[j]!=0)||(j==(length-1))))
            {
                if(site[j]==0)
                    temp_end = j;
                else
                    temp_end = j-1;

                temp_mode = mode;
                if(temp_start == 0)
                    start_site = start;
                else
                {
                    start_site = site[temp_start-1];
                    if(temp_mode == 1)
                        temp_mode = 2;
                }

                if(j == length-1)
                    end_site = end;
                else
                {
                    if(site[temp_end+1]>end)
                        end_site = end;
                    else
                        end_site = site[temp_end+1];
                    if(temp_mode == 0)
                        temp_mode = 2;
                }

                length_t = temp_end-temp_start+1;
                char stext[MAX_READ_LENGTH];
                strncpy(stext,&seq[temp_start],length_t);
                stext[length_t] = '\0';
                for(k = 0; k<length_t; k++)
                    stext[k] = toupper(stext[k]);

                if(end_site<=start_site+1)//前后尾端位置偏差
                {
                    if((mode==1)&&(temp_end==(length-1)))
                    {
                        if(end_site<start_site+1)
                        {
                            i = 1;
                            while(((temp_start-i)>=0)&&(site[temp_start-i] > end_site))
                            {
                                site[temp_start-i] = end_site-1;
                                i++;
                            }
                        }
                        for(i = temp_start;i<=temp_end;i++)
                            site[i] = end_site-1;
                    }
                    else
                    {
                        if(end_site<start_site+1)
                        {
                            i = 1;
                            while(((temp_end+i)<length)&&(site[temp_end+i] < start_site))
                            {
                                site[temp_end+i] = start_site;
                                i++;
                            }
                        }
                        for(i = temp_start;i<=temp_end;i++)
                            site[i] = start_site;
                    }

                    flag = 0;
                    continue;
                }
                if (find8(chr,stext,temp_mode,start_site,end_site,hash,hash_num,seed,&(site[temp_start]),1,chr_order,DP)==1) {return 1;}
                flag = 0;
            }
        }
    }
    {return 0;}
};
int find8_l(char *chr,char *seq,int mode,unsigned int start,unsigned int end,
               unsigned int **hash,int *hash_num,struct seed_t *seed,unsigned int *site,int align_flag,int chr_order,struct Splice_DP_t *DP)
{
    int i,j,k;
    int length = strlen(seq);

    int length_t;

    char f_seq[KMER_FRONT+1];
    uint32_t f_base = 0;
    uint32_t fmax = pow(4,KMER_FRONT);
    int hash_site = 0;

    int seed_num = 0;

    int temp_start = 0;
    int temp_end = 0;
    int temp_mode = 0;
    unsigned int start_site = 0;
    unsigned int end_site = 0;
    int flag = 0;

    int step = 0;

    unsigned int chr_start = opt->chr->list[chr_order].start_site;
    unsigned int START = start+chr_start;
    unsigned int END = end+chr_start;

    {
    {
        step = 4;//KMER_FRONT/2;
        i = 0;
        while ((i+KMER_FRONT)<=length)
        {
            strncpy(f_seq,&seq[i],KMER_FRONT);
            f_seq[KMER_FRONT] = '\0';
            f_base = char2base_16(f_seq,KMER_FRONT);

            if(f_base>=fmax) {return -1;}

            hash_site = hash_c_find(hash,hash_num,f_base,START);
            if((hash_site+AREA_SEED_CAND_NUM<hash_num[f_base])&&(hash[f_base][hash_site+AREA_SEED_CAND_NUM]>=START)&&(hash[f_base][hash_site+AREA_SEED_CAND_NUM]<=END)) {i+=step;continue;}
            while ((hash_site < hash_num[f_base])&&(hash[f_base][hash_site]<=END))
            {
                if((hash[f_base][hash_site]>=START)&&(hash[f_base][hash_site]<=END))
                {
                    seed[seed_num].pos = hash[f_base][hash_site]-chr_start;
                    seed[seed_num].start = i;
                    seed[seed_num].length = KMER_FRONT;
                    seed[seed_num].abs = seed[seed_num].pos-seed[seed_num].start;
                    seed[seed_num].last = -1;
                    seed[seed_num].score = 0;
                    seed_num++;
                    //insert_seed(seed,&seed_num,SEED_BUF_LENGTH,seed_r);
                }
                hash_site++;
            }
            i+=step;
        }
    }


        int s = 0;
        int l = 0;
        //seed前后延伸
        for(i = 0;i<seed_num;i++)
        {
            k = 0;
            s = seed[i].start;
            for(j = 1;j<=s;j++)
            {
                if(nst_nt4_table[(int)chr[seed[i].pos-j]]!=nst_nt4_table[(int)seq[s-j]])
                    break;
                k = j;
            }
            seed[i].pos-=k;
            seed[i].start-=k;
            seed[i].length+=k;

            s = seed[i].start;
            l = seed[i].length;
            k = 0;
            for(j = 0;j<length-s-l;j++)
            {
                if(nst_nt4_table[(int)chr[seed[i].pos+j+l]]!=nst_nt4_table[(int)seq[s+l+j]])
                    break;
                k = j+1;
            }
            seed[i].length+=k;
        }
        qsort(seed,seed_num,sizeof(struct seed_t),seed_cmp);

        int seed_start = 0;
        for(i = 0;i<seed_num;i++)
        {
            while((seed[i].start >seed[seed_start].start)&&(seed_start<seed_num)) seed_start++;
            for(j = seed_start;j<i;j++)
            {
                if(seed[j].length == 0)continue;
                if((seed[i].start == seed[j].start)&&(seed[i].pos == seed[j].pos))
                {seed[i].start = 0;seed[i].pos = 0;seed[i].length = 0;seed[i].score = 0;seed[i].last = -1;break;}
            }
        }
        qsort(seed,seed_num,sizeof(struct seed_t),seed_cmp);
        int init_num = seed_num;
        extend_seed_tail(chr,chr_order,seq,mode,start,end,seed,&seed_num,init_num);

        int max_o = -1;
        int max_s = 0;

        if(mode!=1) qsort(seed,seed_num,sizeof(struct seed_t),seed_cmp);
        else qsort(seed,seed_num,sizeof(struct seed_t),seed_cmp_r);

        seed_start = 0;
        for(i = 0;i<seed_num;i++)
        {
            while((seed[i].start >seed[seed_start].start)&&(seed_start<seed_num)) seed_start++;
            for(j = seed_start;j<i;j++)
            {
                if(seed[j].length == 0)continue;
                if((seed[i].start == seed[j].start)&&(seed[i].pos == seed[j].pos))
                {seed[i].start = 0;seed[i].pos = 0;seed[i].length = 0;seed[i].score = 0;seed[i].last = -1;break;}
            }
        }

        for(i = 0; i<seed_num; i++)
        {
            if(seed[i].length == 0){seed[i].score = -1000;continue;}

            seed[i].score = seed[i].length*(opt->match);
            if(mode!=1)
            {
                if((seed[i].abs>=start+1)&&(seed[i].abs-start>opt->change_length)) seed[i].score-=opt->splice;
                //else if ((seed[i].abs<start+1)&&(start+1-seed[i].abs>=15)) {seed[i].score = -1000;continue;}
                else if ((seed[i].abs>=start+1)&&(seed[i].abs-start-1<=opt->change_length)) seed[i].score-=4+(seed[i].abs-start-1-1)*opt->gap;
                else if ((seed[i].abs<start+1)) seed[i].score-=4+(start+1-seed[i].abs)*opt->match+(start+1-seed[i].abs-1)*opt->gap;
            }

            max_o = -1;
            max_s = 0;
            for(j = 0;j<i;j++)
            {
                if((seed[i].start > seed[j].start)&&((seed[i].start+seed[i].length)>(seed[j].start+seed[j].length)&&(seed[i].pos > seed[j].pos)))
                {
                    if(seed[i].start<=seed[j].start+seed[j].length)
                        max_s = seed[j].score + (seed[i].start+seed[i].length-seed[j].start-seed[j].length)*opt->match;
                    else
                        max_s = seed[j].score + seed[i].length*opt->match;

                    if(seed[i].abs<seed[j].abs)// max_s--;//max_s -=(seed[j].abs-seed[i].abs)+1;
                    max_s -=(seed[j].abs-seed[i].abs)*opt->match+4+(seed[j].abs-seed[i].abs-1)*opt->gap;
                    //{
                        //l = seed[j].abs-seed[i].abs;
                        //if(l>15) continue;
                       // else max_s -=(seed[j].abs-seed[i].abs)*opt->match+4+(seed[j].abs-seed[i].abs-1)*opt->gap;
                    //}
                    else if((seed[i].abs>seed[j].abs)&&(seed[i].abs-seed[j].abs<opt->change_length))
                        max_s -=4+(seed[i].abs-seed[j].abs-1)*opt->gap;
                    else if((seed[i].abs>seed[j].abs)&&(seed[i].abs-seed[j].abs>opt->change_length))
                        max_s -=opt->splice;

                    if(seed[i].score < max_s)
                    {
                        max_o = j;
                        seed[i].score = max_s;
                    }
                }
            }
            if(max_o != -1) seed[i].last = max_o;
            //if(max_o == -1)
            //{
                //seed[i].score = seed[i].length*(opt->match);
                //if((mode!=1)&&(seed[i].abs-start>opt->change_length)) seed[i].score-=opt->splice;
            //}
            //else
                //seed[i].last = max_o;

        }
        //回溯
        max_o = -1;
        max_s = -1000;
        int best_num = 0;
        for(i = seed_num-1; i>=0; i--)
        {
            if(mode!=0)
            {
                if((end-length>seed[i].abs)&&(end-length-seed[i].abs>opt->change_length)) seed[i].score-=opt->splice;
                //else if((end-length<seed[i].abs)&&(seed[i].abs-(end-length)>15)) continue;

                else if((end-length>seed[i].abs)&&(end-length-seed[i].abs<=opt->change_length)) seed[i].score-=4+(end-length-seed[i].abs-1)*opt->gap;
                else if((end-length<seed[i].abs)) seed[i].score-=4+(seed[i].abs-(end-length))*opt->match+(seed[i].abs-(end-length)-1)*opt->gap;
            }
            if(seed[i].score>=max_s)
            {
                max_o = i;
                max_s = seed[i].score;
                best_num = 1;
            }
        }
        if(best_num>1) max_o = -1;
        //if(max_s<0) max_o = -1;
        while(max_o!=-1)
        {
            for(i = 0;i<seed[max_o].length;i++)
                site[seed[max_o].start+i] = seed[max_o].pos + i;
            if(seed[max_o].last!=-1)
                max_o = seed[max_o].last;
            else max_o = -1;
        }
        for(j = 0;j<length;j++)
        {
            if(site[j]>=end) site[j] = end-1;
            else if((site[j]!=0)&&(site[j]<start)) site[j] = start;

            if((j!=0)&&(site[j]!=0)&&(site[j]<site[j-1]))
                site[j]=site[j-1];

            if ((flag==0)&&(site[j]==0))
            {
                temp_start = j;
                flag = 1;
            }
            if ((flag==1)&&((site[j]!=0)||(j==(length-1))))
            {
                if(site[j]==0)
                    temp_end = j;
                else
                    temp_end = j-1;

                temp_mode = mode;
                if(temp_start == 0)
                    start_site = start;
                else
                {
                    start_site = site[temp_start-1];
                    if(temp_mode == 1)
                        temp_mode = 2;
                }

                if(j == length-1)
                    end_site = end;
                else
                {
                    if(site[temp_end+1]>end)
                        end_site = end;
                    else
                        end_site = site[temp_end+1];
                    if(temp_mode == 0)
                        temp_mode = 2;
                }

                length_t = temp_end-temp_start+1;
                char stext[MAX_READ_LENGTH];
                strncpy(stext,&seq[temp_start],length_t);
                stext[length_t] = '\0';
                for(k = 0; k<length_t; k++)
                    stext[k] = toupper(stext[k]);

                if(end_site<=start_site+1)//前后尾端位置偏差
                {
                    if((mode==1)&&(temp_end==(length-1)))
                    {
                        if(end_site<start_site+1)
                        {
                            i = 1;
                            while(((temp_start-i)>=0)&&(site[temp_start-i] > end_site))
                            {
                                site[temp_start-i] = end_site-1;
                                i++;
                            }
                        }
                        for(i = temp_start;i<=temp_end;i++)
                            site[i] = end_site-1;
                    }
                    else
                    {
                        if(end_site<start_site+1)
                        {
                            i = 1;
                            while(((temp_end+i)<length)&&(site[temp_end+i] < start_site))
                            {
                                site[temp_end+i] = start_site;
                                i++;
                            }
                        }
                        for(i = temp_start;i<=temp_end;i++)
                            site[i] = start_site;
                    }

                    flag = 0;
                    continue;
                }
                //if (length_t> 500) {return 1;}
                if (find8(chr,stext,temp_mode,start_site,end_site,hash,hash_num,seed,&(site[temp_start]),0,chr_order,DP)==1) {return 1;}
                flag = 0;
            }
        }
    }
    {return 0;}
};
int area_align(char *chr,char *seq,int mode,unsigned int start,unsigned int end,
               unsigned int **hash,int *hash_num,struct seed_t *seed,unsigned int *site,int chr_order,struct NW_list *NW,struct Splice_DP_t *DP)//mode 0前为准 1后为准 2前后
{
    if(end<=start) return 0;

    int length = strlen(seq);
    int i = 0;//,k = 0;

    int flag = 0;

    int score = 0,lscore = -10000;
    //int length_r ;

    unsigned int *lsite = NW->site;

    struct cigar_t *t_cigar = (struct cigar_t *)malloc((MAX_CIGAR_BUF+2)*sizeof(struct cigar_t));
    int t_num;

    unsigned int start_pos = 0;

    /*if((end-start-1)<=2*max(length*(1.1),length+20))
    {
        Area_DP(chr,start+1,end-1,seq,length,t_cigar,&t_num,DP);
        score = cigar_score(t_cigar,t_num);
        cigar2site(start+1,site,t_cigar,t_num);
        return score;
    }*/

    if((mode==0)&&(length<200))
    {
        Tail_DP(chr,start+1,end,1,seq,length,t_cigar,&t_num,DP,&start_pos,1);
        lscore = cigar_score(t_cigar,t_num);
        cigar2site(start+1,lsite,t_cigar,t_num);
    }
    else if((mode==1)&&(length<200))
    {
        Tail_DP(chr,end-1,start,0,seq,length,t_cigar,&t_num,DP,&start_pos,1);
        lscore = cigar_score(t_cigar,t_num);
        cigar2site2(end-1,lsite,length,t_cigar,t_num);
    }
    else if((mode==2)&&(length<50))
    {
        if((end-start-1)<=2*max(length*1.1,length+20))
            Area_DP(chr,start+1,end-1,seq,length,t_cigar,&t_num,DP);
        else
            Splice_DP(chr,start+1,end-1,seq,length,t_cigar,&t_num,DP);
        lscore = cigar_score(t_cigar,t_num);
        cigar2site(start+1,lsite,t_cigar,t_num);
    }
    free(t_cigar);

    if (length>500)
    {
        if(find8_l(chr,seq,mode,start,end,hash,hash_num,seed,site,0,chr_order,DP)==1) {return -1000;}
    }
    else if(find8(chr,seq,mode,start,end,hash,hash_num,seed,site,0,chr_order,DP)==1) {return -1000;}

    int out_flag = 0;
    score = 0;
    flag = 1;
    if(mode!=1)
    {
        if((site[0]-start>1)&&(site[0]-start<opt->change_length))score-=opt->gap;
        else if (site[0]-start>opt->change_length) score-=opt->splice;

        if(start==site[0]) {score-=4;flag = 0;}
        else
        {
            if(nst_nt4_table[(int)chr[site[0]]]==nst_nt4_table[(int)seq[0]]) score+=opt->match;
            else score-=opt->miss;
        }
    }
    else
    {
        if(nst_nt4_table[(int)chr[site[0]]]==nst_nt4_table[(int)seq[0]]) score+=opt->match;
        else score-=opt->miss;
    }
    if (site[0]==0) {score=0;out_flag = 1;}
    for(i = 1;i<length;i++)
    {
        if (site[i]==0) {score=0;out_flag = 1;break;}

        if(site[i]-site[i-1]==0)
        {
            if(flag==1) {score-=4;flag = 0;}
            else score-=opt->gap;
        }
        else
        {
            if(nst_nt4_table[(int)chr[site[i]]]==nst_nt4_table[(int)seq[i]]) score+=opt->match;
            else score-=opt->miss;
            flag = 1;
        }

        if((site[i]-site[i-1]>1)&&(site[i]-site[i-1]<=opt->change_length))
        {
            flag = 1;
            score-=4+(site[i]-site[i-1]-2)*opt->gap;
        }
        else if(site[i]-site[i-1]>opt->change_length)
        {
            flag = 1;
            score-=opt->splice;
        }
    }
    if(mode!=0)//if((mode!=0)&&(out_flag==0))
    {
        if((end-site[length-1]>1)&&(end-site[length-1]<opt->change_length))score-=4+(end-site[length-1]-1)*opt->gap;
        else if (end-site[length-1]>opt->change_length) score-=opt->splice;
        else if((end==site[length-1])&&(flag)) score-=4;
        else if((end==site[length-1])&&(!flag)) score-=opt->gap;
    }
    if(lscore>=score)//if((lscore>=score)||(out_flag==1))//((mode!=2)&&((lscore>score)||(out_flag==1)))
    {
        memcpy(site,lsite,sizeof(unsigned int)*length);
        score = lscore;
    }
    if((mode==2)&&(out_flag==1)) return -1000;
    return score;
}
int site2cigar(char *chr,char *seq,unsigned int *site,int length,unsigned int start,unsigned int end,struct cigar_t cigar[MAX_CIGAR_BUF],int *cigar_num,int max,int mode)
{
    int change_length;
    int snp_num = 0;
    int i;

    if (site[0]==0)
    {
        cigar[(*cigar_num)].c = 'S';
        cigar[(*cigar_num)].l = 1;
    }
    else
    {
    if(mode!=1)
    {
        if((site[0]-start>1)&&(site[0]-start<=opt->change_length+1))
        {
            cigar[(*cigar_num)].c = 'D';
            cigar[(*cigar_num)].l = site[0]-start-1;
            (*cigar_num)++;
            if((*cigar_num)>=max) {(*cigar_num)=0;return -1000;}
            snp_num++;
        }
        else if (site[0]-start>opt->change_length+1)
        {
            cigar[(*cigar_num)].c = 'N';
            cigar[(*cigar_num)].l = site[0]-start-1;
            (*cigar_num)++;
            if((*cigar_num)>=max) {(*cigar_num)=0;return -1000;}
        }

        if(start==site[0])
        {
            cigar[(*cigar_num)].c = 'I';
            cigar[(*cigar_num)].l = 1;
            snp_num++;
        }
        else
        {
            if(nst_nt4_table[(int)chr[site[0]]]==nst_nt4_table[(int)seq[0]])
            {
                cigar[(*cigar_num)].c = 'M';
                cigar[(*cigar_num)].l = 1;
            }
            else
            {
                cigar[(*cigar_num)].c = 'X';
                cigar[(*cigar_num)].l = 1;
                snp_num++;
            }
        }
    }
    else
    {
        if(nst_nt4_table[(int)chr[site[0]]]==nst_nt4_table[(int)seq[0]])
        {
            cigar[(*cigar_num)].c = 'M';
            cigar[(*cigar_num)].l = 1;
        }
        else
        {
            cigar[(*cigar_num)].c = 'X';
            cigar[(*cigar_num)].l = 1;
            snp_num++;
        }
    }

    }

    for(i = 1;i<length;i++)
    {
        if (site[i]==0)
        {
            if(cigar[(*cigar_num)].c!='S')
            {
                (*cigar_num)++;
                if((*cigar_num)>=max) {(*cigar_num)=0;return -1000;}
                cigar[(*cigar_num)].c = 'S';
                cigar[(*cigar_num)].l = 1;
            }
            else cigar[(*cigar_num)].l++;
            continue;
        }

        change_length = site[i] - site[i-1]-1;
        if((cigar[(*cigar_num)].c != 'S')&&(change_length > 0)&&(change_length<=opt->change_length))
        {
            (*cigar_num)++;
            if((*cigar_num)>=max) {(*cigar_num)=0;return -1000;}
            cigar[(*cigar_num)].c = 'D';
            cigar[(*cigar_num)].l = change_length;
            snp_num++;
        }
        else if((cigar[(*cigar_num)].c != 'S')&&(change_length>opt->change_length))
        {
            (*cigar_num)++;
            cigar[(*cigar_num)].c = 'N';
            cigar[(*cigar_num)].l = change_length;
            snp_num++;
        }

        if(site[i] == site[i-1])
        {
            if(cigar[(*cigar_num)].c!='I')
            {
                (*cigar_num)++;
                if((*cigar_num)>=max) {(*cigar_num)=0;return -1000;}
                cigar[(*cigar_num)].c = 'I';
                cigar[(*cigar_num)].l = 1;
                snp_num++;
            }
            else cigar[(*cigar_num)].l++;
        }
        else
        {
            if(nst_nt4_table[(int)seq[i]]==nst_nt4_table[(int)chr[site[i]]])
            {
                if(cigar[(*cigar_num)].c!='M')
                {
                    (*cigar_num)++;
                    if((*cigar_num)>=max) {(*cigar_num)=0;return -1000;}

                    cigar[(*cigar_num)].c = 'M';
                    cigar[(*cigar_num)].l = 1;
                }
                else cigar[(*cigar_num)].l++;
            }
            else
            {
                if(cigar[(*cigar_num)].c!='X')
                {
                    (*cigar_num)++;
                    if((*cigar_num)>=max) {(*cigar_num)=0;return -1000;}

                    cigar[(*cigar_num)].c = 'X';
                    cigar[(*cigar_num)].l = 1;
                }
                else cigar[(*cigar_num)].l++;
                snp_num++;
            }
        }
    }
    (*cigar_num)++;
    if((*cigar_num)>=max) {(*cigar_num)=0;return -1000;}

    if(mode!=0)
    {
        if((cigar[(*cigar_num)-1].c != 'S')&&(end-site[length-1]>1)&&(end-site[length-1]<=opt->change_length+1))
        {
            cigar[(*cigar_num)].c = 'D';
            cigar[(*cigar_num)].l = end-site[length-1]-1;
            (*cigar_num)++;
            if((*cigar_num)>=max) {(*cigar_num)=0;return -1000;}
            snp_num++;
        }
        else if((cigar[(*cigar_num)-1].c != 'S')&&(end-site[length-1]>opt->change_length+1))
        {
            cigar[(*cigar_num)].c = 'N';
            cigar[(*cigar_num)].l = end-site[length-1]-1;
            (*cigar_num)++;
            if((*cigar_num)>=max) {(*cigar_num)=0;return -1000;}
            snp_num++;
        }
    }
    return snp_num;
}
void write_cigar(struct cigar_t *cigar,int cigar_num,char *cigarS)
{
    int i;
    char temp[20];
    cigarS[0] = '\0';
    int M = 0;
    for(i = 0;i<cigar_num;i++)
    {
        if((cigar[i].c=='M')||(cigar[i].c=='X')) M+=cigar[i].l;
        else
        {
            if(M!=0)
            {
                sprintf(temp,"%dM",M);
                strcat(cigarS,temp);
                M=0;
            }
            if((cigar[i].l<=opt->change_length)&&(cigar[i].c=='N')) cigar[i].c = 'D';
            sprintf(temp,"%d%c",cigar[i].l,cigar[i].c);
            strcat(cigarS,temp);
        }
    }
    if(M!=0)
    {
        sprintf(temp,"%dM",M);
        strcat(cigarS,temp);
    }
}
void generate_cigar(struct read_t *read,struct seed_t *seed,int *seed_order,int seed_order_n,struct seed_t *seed_a,struct result_t *result,int *result_num,int read_order,FILE *out,FILE *un,
                    struct NW_list *NW1,struct NW_list *NW2,struct NW_list *NW3,struct Splice_DP_t *DP,struct splice_list *SP,int *SP_num)
{
    int i = 0,j = 0,x = 0,r = 0;
    uint64_t ref_pos = 0,start_pos = 0;
    int read_pos;

    int chr_order = 0;
    unsigned int chr_start = 0;
    char strand = 0;

    struct cigar_t *cigar = (struct cigar_t *)malloc(MAX_CIGAR_BUF*sizeof(struct cigar_t));
    int cigar_num = 0;

    struct cigar_t *t_cigar = (struct cigar_t *)malloc(MAX_CIGAR_BUF*sizeof(struct cigar_t));
    int t_num;

    //int snp_num;

    int break_flag = 0;
    struct seed_t temp_seed;

    unsigned int *site = NW1->site;
    char *text = NW1->text;

    unsigned int *tsite = NW2->site;
    char *ttext = NW2->text;

    int t;
    int l;

    int front;

    if(seed[seed_order[0]].pos>=opt->l_pac)
    {
        strand = 1;

        int temp_n = 0;
        for(i = 0;i<seed_order_n;i++)
        {
            seed[seed_order[i]].pos = (opt->l_pac<<1)-seed[seed_order[i]].pos-seed[seed_order[i]].length;
            seed[seed_order[i]].start = read->length-seed[seed_order[i]].start-seed[seed_order[i]].length;
            seed[seed_order[i]].abs = seed[seed_order[i]].pos - seed[seed_order[i]].start;
        }
        for(i = 0;i<seed_order_n/2;i++)
        {
            temp_n = seed_order[i];
            seed_order[i] = seed_order[seed_order_n-1-i];
            seed_order[seed_order_n-1-i] = temp_n;
        }
    }

    for(j = 0;j<opt->chr->total;j++)
    {
        if((seed[seed_order[0]].pos >= opt->chr->list[j].start_site)
        &&(seed[seed_order[0]].pos < opt->chr->list[j].start_site+opt->chr->list[j].length))
        {
            chr_order = j;
            chr_start = opt->chr->list[chr_order].start_site;
            break;
        }
    }
    if ((chr_order>=opt->chr->total)||(chr_order<0))
        break_flag = 1;

    cigar_num = 0;

    ref_pos = seed[seed_order[seed_order_n-1]].pos;
    start_pos = ref_pos-chr_start;
    read_pos = 0;
    memcpy(&temp_seed,&(seed[seed_order[seed_order_n-1]]),sizeof(struct seed_t));
    if(seed[seed_order[seed_order_n-1]].start!=0)
    {
        /*cigar[cigar_num].c = 'S';
        cigar[cigar_num].l = seed[seed_order[seed_order_n-1]].start;
        cigar_num++;*/

        if(seed[seed_order[seed_order_n-1]].start<=20)
        {
        l = min(5,seed[seed_order[seed_order_n-1]].length);
        for(j = 0;j<seed[seed_order[seed_order_n-1]].start+l;j++)
        {
            if(strand==0) ttext[j] = read->seq[j];
            else ttext[j] = read->rseq[j];
        }
        ttext[j] = '\0';
        memset(tsite, 0, (seed[seed_order[seed_order_n-1]].start+l)*sizeof(unsigned int)/sizeof(char));
        t=area_align(opt->chr->list[chr_order].seq,ttext,1,((ref_pos>chr_start+AREA)?(ref_pos+l-chr_start-AREA):0),ref_pos+l-chr_start,opt->c_hash,opt->c_num,seed_a,&(tsite[0]),chr_order,NW3,DP);
        }
        else t = -10000;

        //if(seed[seed_order[seed_order_n-1]].start>20)
        {
        for(j = 0;j<seed[seed_order[seed_order_n-1]].start;j++)
        {
            if(strand==0) text[j] = read->seq[j];
            else text[j] = read->rseq[j];
        }
        text[j] = '\0';
        memset(site, 0, seed[seed_order[seed_order_n-1]].start*sizeof(unsigned int)/sizeof(char));
        r=area_align(opt->chr->list[chr_order].seq,text,1,((ref_pos>chr_start+AREA)?(ref_pos-chr_start-AREA):0),ref_pos-chr_start,opt->c_hash,opt->c_num,seed_a,&(site[0]),chr_order,NW3,DP);
        }
        //else r = -1000;

        if(r+l*opt->match<t)
        {
            j = 0;
            while(j<seed[seed_order[seed_order_n-1]].start)
            {
                if(tsite[j]!=0) {start_pos = tsite[j];break;}
                j++;
            }
            t_num = 0;
            l = min(5,seed[seed_order[seed_order_n-1]].length);
            site2cigar(opt->chr->list[chr_order].seq,ttext,&(tsite[0]),seed[seed_order[seed_order_n-1]].start+l,((ref_pos>chr_start+AREA)?(ref_pos+l-chr_start-AREA):0),ref_pos+l-chr_start,t_cigar,&t_num,MAX_CIGAR_BUF,1);
            {
                if(cigar_num+t_num>MAX_CIGAR_BUF){break_flag = 1;}
                else
                     {memcpy(cigar,t_cigar,t_num*sizeof(struct cigar_t));cigar_num+=t_num;}
                cigar[cigar_num].l = seed[seed_order[seed_order_n-1]].length-l;
                cigar[cigar_num].c = 'M';
                if(cigar[cigar_num].l==0) cigar_num--;
            }
        }
        else
        {
        if((r==-1000))//||(seed[seed_order[seed_order_n-1]].start<20))
        {
            cigar[cigar_num].c = 'S';
            cigar[cigar_num].l = seed[seed_order[seed_order_n-1]].start;
            cigar_num++;
            cigar[cigar_num].l = seed[seed_order[seed_order_n-1]].length;
            cigar[cigar_num].c = 'M';
        }
        else
        {
            j = 0;
            while(j<seed[seed_order[seed_order_n-1]].start)
            {
                if(site[j]!=0) {start_pos = site[j];break;}
                j++;
            }
            t_num = 0;
            site2cigar(opt->chr->list[chr_order].seq,text,&(site[0]),seed[seed_order[seed_order_n-1]].start,((ref_pos>chr_start+AREA)?(ref_pos-chr_start-AREA):0),ref_pos-chr_start,t_cigar,&t_num,MAX_CIGAR_BUF,1);
            {
                if(cigar_num+t_num>MAX_CIGAR_BUF){break_flag = 1;}
                else
                     {memcpy(cigar,t_cigar,t_num*sizeof(struct cigar_t));cigar_num+=t_num;}
            }
            cigar[cigar_num].l = seed[seed_order[seed_order_n-1]].length;
            cigar[cigar_num].c = 'M';
        }
        }
    }
    else
    {
        cigar[cigar_num].l = seed[seed_order[seed_order_n-1]].length;
        cigar[cigar_num].c = 'M';
    }

    ref_pos = seed[seed_order[seed_order_n-1]].pos+seed[seed_order[seed_order_n-1]].length;
    read_pos = seed[seed_order[seed_order_n-1]].start+seed[seed_order[seed_order_n-1]].length;

    //check_cigar2(cigar,cigar_num+1,read_pos);

    front = seed_order_n-1;
    for (i = seed_order_n-2;i>=0;i--)
    {
        if(seed[seed_order[i]].abs==seed[seed_order[front]].abs)
        {
            if(seed[seed_order[i]].start<=read_pos)//基本不可能
            {
                cigar[cigar_num].l += seed[seed_order[i]].length-(read_pos-seed[seed_order[i]].start);

                ref_pos+=seed[seed_order[i]].length-(read_pos-seed[seed_order[i]].start);
                read_pos+=seed[seed_order[i]].length-(read_pos-seed[seed_order[i]].start);
            }
            else //XM *
            {
                x = seed[seed_order[i]].start-read_pos;
                for(j = 0;j<x;j++)
                {
                    if(((strand==0)&&(nst_nt4_table[(int)read->seq[read_pos+j]]==nst_nt4_table[(int)opt->chr->list[chr_order].seq[ref_pos+j-chr_start]]))
                        ||((strand==1)&&(nst_nt4_table[(int)read->rseq[read_pos+j]]==nst_nt4_table[(int)opt->chr->list[chr_order].seq[ref_pos+j-chr_start]])))
                    {
                        if(cigar[cigar_num].c=='M')
                            cigar[cigar_num].l++;
                        else
                        {
                            cigar_num++;
                            if(cigar_num>MAX_CIGAR_BUF){break_flag = 1;break;}
                            cigar[cigar_num].c='M';
                            cigar[cigar_num].l = 1;
                        }
                    }
                    else
                    {
                        if(cigar[cigar_num].c=='X')
                            cigar[cigar_num].l++;
                        else
                        {
                            cigar_num++;
                            if(cigar_num>MAX_CIGAR_BUF){break_flag = 1;break;}
                            cigar[cigar_num].c='X';
                            cigar[cigar_num].l = 1;
                        }
                    }
                }
                cigar_num++;
                if(cigar_num>MAX_CIGAR_BUF){break_flag = 1;break;}
                cigar[cigar_num].l = seed[seed_order[i]].length;
                cigar[cigar_num].c = 'M';

                ref_pos = seed[seed_order[i]].pos+seed[seed_order[i]].length;
                read_pos = seed[seed_order[i]].start+seed[seed_order[i]].length;
            }
        }
        else if(seed[seed_order[i]].abs<=seed[seed_order[front]].abs)
        {
            cigar_num++;
            if(cigar_num>MAX_CIGAR_BUF){break_flag = 1;break;}

            memcpy(&temp_seed,&(seed[seed_order[i]]),sizeof(struct seed_t));
            if(seed[seed_order[i]].start>read_pos) //IX *
            {
                if(seed[seed_order[i]].pos<=ref_pos)
                {
                    cigar[cigar_num].l = seed[seed_order[front]].abs-seed[seed_order[i]].abs;
                    cigar[cigar_num].c = 'I';

                    read_pos+=cigar[cigar_num].l;
                    temp_seed.pos +=read_pos-seed[seed_order[i]].start;
                    temp_seed.length -=read_pos-seed[seed_order[i]].start;
                    temp_seed.start +=read_pos-seed[seed_order[i]].start;

                    if((temp_seed.length<=0))//||(cigar[cigar_num].l>opt->change_length))
                        {read_pos-=cigar[cigar_num].l;cigar_num--;continue;}

                    cigar_num++;
                    if(cigar_num>MAX_CIGAR_BUF){break_flag = 1;break;}
                }
                else
                {
                    for(j = 0;j<seed[seed_order[i]].start-read_pos;j++)
                    {
                        if(strand==0) text[j] = read->seq[j+read_pos];
                        else text[j] = read->rseq[j+read_pos];
                    }
                    text[j] = '\0';

                    //if(seed[seed_order[i]].pos-ref_pos <= MAX_DP_AREA)
                    //{
                        //Area_Splice_DP(opt->chr->list[chr_order].seq,ref_pos-chr_start,seed[seed_order[i]].pos-chr_start-1,text,j,t_cigar,&t_num,DP);
                        //{memcpy(&(cigar[cigar_num]),t_cigar,t_num*sizeof(struct cigar_t));cigar_num+=t_num;}
                    //}
                    //else
                    {

                    memset(site, 0, (seed[seed_order[i]].start-read_pos)*sizeof(unsigned int)/sizeof(char));
                    r=area_align(opt->chr->list[chr_order].seq,text,2,ref_pos-1-chr_start,seed[seed_order[i]].pos-chr_start,opt->c_hash,opt->c_num,seed_a,&(site[0]),chr_order,NW3,DP);
                    if(r==-1000)
                    {
                        cigar[cigar_num].c = 'S';
                        cigar[cigar_num].l = seed[seed_order[i]].start-read_pos;
                        cigar_num++;

                        cigar[cigar_num].c = 'N';
                        cigar[cigar_num].l = seed[seed_order[i]].pos-ref_pos;
                        cigar_num++;
                    }
                    else
                    {
                        t_num = 0;
                        site2cigar(opt->chr->list[chr_order].seq,text,&(site[0]),seed[seed_order[i]].start-read_pos,ref_pos-1-chr_start,seed[seed_order[i]].pos-chr_start,t_cigar,&t_num,MAX_CIGAR_BUF,2);
                        if(cigar_num+t_num>MAX_CIGAR_BUF){break_flag = 1;break;}
                        else
                        {
                            if(cigar[cigar_num-1].c==t_cigar[0].c)
                            {
                                cigar[cigar_num-1].l+=t_cigar[0].l;
                                memcpy(cigar+cigar_num,t_cigar+1,(t_num-1)*sizeof(struct cigar_t));
                                cigar_num+=t_num-1;
                            }
                            else
                            {memcpy(cigar+cigar_num,t_cigar,t_num*sizeof(struct cigar_t));cigar_num+=t_num;}
                        }
                    }

                    }
                }
            }
            else
            {
                cigar[cigar_num].l = seed[seed_order[front]].abs-seed[seed_order[i]].abs;
                cigar[cigar_num].c = 'I';

                read_pos+=cigar[cigar_num].l;
                temp_seed.pos +=read_pos-seed[seed_order[i]].start;
                temp_seed.length -=read_pos-seed[seed_order[i]].start;
                temp_seed.start +=read_pos-seed[seed_order[i]].start;

                if((temp_seed.length<=0))//||(cigar[cigar_num].l>opt->change_length))
                    {read_pos-=cigar[cigar_num].l;cigar_num--;continue;}

                cigar_num++;
                if(cigar_num>MAX_CIGAR_BUF){break_flag = 1;break;}
            }
            if(temp_seed.length>0)
            {
                if(cigar[cigar_num-1].c=='M')
            {
                cigar_num--;
                cigar[cigar_num].l+=temp_seed.length;
            }
            else
            {
                cigar[cigar_num].l = temp_seed.length;
                cigar[cigar_num].c = 'M';
            }
                //cigar[cigar_num].l = temp_seed.length;
                //cigar[cigar_num].c = 'M';

                ref_pos = seed[seed_order[i]].pos+seed[seed_order[i]].length;
                read_pos = seed[seed_order[i]].start+seed[seed_order[i]].length;
            }
            else {
                if(i==0)
                {
                    ref_pos -= cigar[cigar_num-1].l;
                    cigar_num--;
                }
                cigar_num--;
                continue;
            }
            //check_cigar2(cigar,cigar_num+1,read_pos);
        }
        else//DN
        {
            cigar_num++;
            if(cigar_num>=MAX_CIGAR_BUF){break_flag = 1;break;}

            memcpy(&temp_seed,&(seed[seed_order[i]]),sizeof(struct seed_t));

            if(read_pos<seed[seed_order[i]].start)
            {
                for(j = 0;j<seed[seed_order[i]].start-read_pos;j++)
                {
                        if(strand==0) text[j] = read->seq[j+read_pos];
                        else text[j] = read->rseq[j+read_pos];
                }
                text[j] = '\0';
                //if(seed[seed_order[i]].pos-ref_pos <= MAX_DP_AREA)
                //{
                   // Area_Splice_DP(opt->chr->list[chr_order].seq,ref_pos-chr_start,seed[seed_order[i]].pos-chr_start-1,text,j,t_cigar,&t_num,DP);
                    //{memcpy(&(cigar[cigar_num]),t_cigar,t_num*sizeof(struct cigar_t));cigar_num+=t_num;}
                //}
                //else
                {
                memset(site, 0, (seed[seed_order[i]].start-read_pos)*sizeof(unsigned int)/sizeof(char));
                r=area_align(opt->chr->list[chr_order].seq,text,2,ref_pos-1-chr_start,seed[seed_order[i]].pos-chr_start,opt->c_hash,opt->c_num,seed_a,&(site[0]),chr_order,NW3,DP);
                if(r==-1000)
                {
                        cigar[cigar_num].c = 'S';
                        cigar[cigar_num].l = seed[seed_order[i]].start-read_pos;
                        cigar_num++;

                        cigar[cigar_num].c = 'N';
                        cigar[cigar_num].l = seed[seed_order[i]].pos-ref_pos;
                        cigar_num++;
                }
                else
                {
                    t_num = 0;
                    site2cigar(opt->chr->list[chr_order].seq,text,&(site[0]),seed[seed_order[i]].start-read_pos,ref_pos-1-chr_start,seed[seed_order[i]].pos-chr_start,t_cigar,&t_num,MAX_CIGAR_BUF,2);
                    if(cigar_num+t_num>MAX_CIGAR_BUF){break_flag = 1;break;}
                    else
                    {
                        if(cigar[cigar_num-1].c==t_cigar[0].c)
                        {
                            cigar[cigar_num-1].l+=t_cigar[0].l;
                            memcpy(cigar+cigar_num,t_cigar+1,(t_num-1)*sizeof(struct cigar_t));
                            cigar_num+=t_num-1;
                        }
                        else
                            {memcpy(cigar+cigar_num,t_cigar,t_num*sizeof(struct cigar_t));cigar_num+=t_num;}
                    }
                }
                }
            }
            else if(read_pos>=seed[seed_order[i]].start)
            {
                temp_seed.pos+=read_pos-seed[seed_order[i]].start;
                temp_seed.length-=read_pos-seed[seed_order[i]].start;
                temp_seed.start+=read_pos-seed[seed_order[i]].start;

                cigar[cigar_num].l = temp_seed.pos-ref_pos;
                if(cigar[cigar_num].l<0)
                {
                    cigar[cigar_num].l = 0-cigar[cigar_num].l;
                    cigar[cigar_num].c = 'I';
                    cigar_num++;
                    if(cigar_num>MAX_CIGAR_BUF){break_flag = 1;break;}
                }
                else if(cigar[cigar_num].l<=opt->change_length)
                {
                    if(cigar[cigar_num-1].c == 'D')
                        {
                            cigar[cigar_num-1].l+=cigar[cigar_num].l;
                            cigar_num--;
                        }
                    cigar[cigar_num].c = 'D';
                    cigar_num++;
                    if(cigar_num>MAX_CIGAR_BUF){break_flag = 1;break;}
                }
                else
                {
                    if(cigar[cigar_num-1].c == 'N')
                        {
                            cigar[cigar_num-1].l+=cigar[cigar_num].l;
                            cigar_num--;
                        }
                    cigar[cigar_num].c = 'N';

                    cigar_num++;
                    if(cigar_num>MAX_CIGAR_BUF){break_flag = 1;break;}
                }
            }
            if(cigar[cigar_num-1].c=='M')
            {
                cigar_num--;
                cigar[cigar_num].l+=temp_seed.length;
            }
            else
            {
                cigar[cigar_num].l = temp_seed.length;
                cigar[cigar_num].c = 'M';
            }

            read_pos=seed[seed_order[i]].start+seed[seed_order[i]].length;
            ref_pos=seed[seed_order[i]].pos+seed[seed_order[i]].length;
            //check_cigar2(cigar,cigar_num+1,read_pos);
        }
        front = i;
    }
    cigar_num++;
    if((read_pos!=read->length)&&(!break_flag))
    {
        if(read_pos>read->length)
            cigar[cigar_num-1].l -= read->length-read_pos;
        else
        {
            if(read->length-read_pos<=20)
            {
            l = min(5,cigar[cigar_num-1].l);

            for(j = 0;j<read->length-read_pos+l;j++)
            {
                if(strand==0) ttext[j] = read->seq[j+read_pos-l];
                else ttext[j] = read->rseq[j+read_pos-l];
            }
            ttext[j] = '\0';
            memset(tsite, 0, (read->length-read_pos+l)*sizeof(unsigned int)/sizeof(char));
            t=area_align(opt->chr->list[chr_order].seq,ttext,0,ref_pos-1-chr_start-l,
                         ((ref_pos+AREA>=opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length)?(opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length-1):(ref_pos-l+AREA))-chr_start,
                         opt->c_hash,opt->c_num,seed_a,&(tsite[0]),chr_order,NW3,DP);
            }
            else t = -10000;

            //if(read->length-read_pos<=20)
            {
            for(j = 0;j<read->length-read_pos;j++)
            {
                if(strand==0) text[j] = read->seq[j+read_pos];
                else text[j] = read->rseq[j+read_pos];
            }
            text[j] = '\0';
            memset(site, 0, (read->length-read_pos)*sizeof(unsigned int)/sizeof(char));
            r=area_align(opt->chr->list[chr_order].seq,text,0,ref_pos-1-chr_start,
                         ((ref_pos+AREA>=opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length)?(opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length-1):(ref_pos+AREA))-chr_start,
                         opt->c_hash,opt->c_num,seed_a,&(site[0]),chr_order,NW3,DP);
            }
            //else r = -1000;

            if(r+l*opt->match<t)
            {
                cigar[cigar_num-1].l-=l;
                t_num = 0;
                site2cigar(opt->chr->list[chr_order].seq,ttext,&(tsite[0]),read->length-read_pos+l,ref_pos-1-chr_start-l,
                         ((ref_pos+AREA>=opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length)?(opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length-1):(ref_pos-l+AREA))-chr_start,
                            t_cigar,&t_num,MAX_CIGAR_BUF,0);
                {
                    if (cigar_num+t_num>MAX_CIGAR_BUF){break_flag = 1;}
                    else
                    {
                        if(t_cigar[0].c==cigar[cigar_num-1].c)
                        {
                            cigar[cigar_num-1].l+=t_cigar[0].l;
                            memcpy(&(cigar[cigar_num]),&(t_cigar[1]),(t_num-1)*sizeof(struct cigar_t));cigar_num+=t_num-1;
                        }
                        else {memcpy(&(cigar[cigar_num]),t_cigar,t_num*sizeof(struct cigar_t));cigar_num+=t_num;}
                    }
                }
            }
            else
            {
            if((r==-1000))//||(read->length-read_pos<20))
            {
                cigar[cigar_num].c = 'S';
                cigar[cigar_num].l = read->length-read_pos;
                cigar_num++;
            }
            else
            {
                t_num = 0;
                site2cigar(opt->chr->list[chr_order].seq,text,&(site[0]),read->length-read_pos,
                ref_pos-1-chr_start,((ref_pos+AREA>=opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length)?(opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length-1):(ref_pos+AREA))-chr_start,
                t_cigar,&t_num,MAX_CIGAR_BUF,0);
                {
                    if (cigar_num+t_num>MAX_CIGAR_BUF){break_flag = 1;}
                    else
                        {memcpy(&(cigar[cigar_num]),t_cigar,t_num*sizeof(struct cigar_t));cigar_num+=t_num;}
                }
            }
            }
        }
    }
    if(break_flag) {cigar_num = 0;}
    else
    {
        char cigarS[MAX_STRING_LENGTH*2];
        check_exon(read-read_order,read->name,(strand==1)?read->rseq:read->seq,read->length,cigar,&cigar_num,&start_pos,chr_order,DP,cigarS,SP,SP_num);
        //check_read_splice(read->seq,read->rseq,cigar,&cigar_num,start_pos,chr_order,strand,DP);

        //check_cigar(read->name,cigarS,read->length);

        int lengthS = strlen(cigarS);

        result[(*result_num)].read_order = read_order;
        result[(*result_num)].cigar = (char *)calloc(lengthS+1,sizeof(char));
        strcpy(result[(*result_num)].cigar,cigarS);
        result[(*result_num)].cigar[lengthS] = '\0';
        result[(*result_num)].chr = chr_order;

        result[(*result_num)].ref_pos  = start_pos;
        result[(*result_num)].strand = strand;
        (*result_num)++;
        if((*result_num)>=opt->result_block)
        {
            pthread_mutex_lock(&OutputLock);
            out_read((read-read_order),result,*result_num,out);
            pthread_mutex_unlock(&OutputLock);
            (*result_num) = 0;
        }
    }
    free(cigar);
    free(t_cigar);
}

