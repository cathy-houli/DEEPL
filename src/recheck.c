#include "main.h"

static pthread_mutex_t ReadLock,OutputLock,UpdateLock;

struct read_inf_t
{
	char *name;
	int length;
	char *seq,*qual;
	int chr;
	uint64_t site;
	char *cigar;
	int flag;
	int NM;
	int score;
	char TAG[100];
};
int find_chr_name(char *chr)
{
    int i = 0;
    for(i = 0;i<opt->chr->total;i++)
    {
        if(strcmp(chr,opt->chr->list[i].name)==0)
            return i;
    }
    return -1;
}
int read_file_sam(struct read_inf_t *read,FILE *file)
{
    int i = 0;
    char *f_line = (char *)calloc(MAX_STRING_LENGTH, sizeof(char));
    char *name = (char *)calloc(MAX_NAME_LENGTH, sizeof(char));
    char *cigar = (char *)calloc(MAX_STRING_LENGTH, sizeof(char));
    char *seq = (char *)calloc(MAX_STRING_LENGTH, sizeof(char));
    char *qual = (char *)calloc(MAX_STRING_LENGTH, sizeof(char));

    while ((i<opt->thread_block)&&(fgets(f_line,MAX_STRING_LENGTH,file)!=NULL))
    {
        if(f_line[0]=='@') continue;
        sscanf(f_line,"%s\t%d\t%d\t%lu\t%*d\t%s\t%*s\t%*u\t%*d\t%s\t%s",name,&read[i].flag,&read[i].chr,&read[i].site,cigar,seq,qual);

        read[i].name = (char *)malloc((strlen(name)+1) *sizeof(char));
        strcpy(read[i].name,name);
        read[i].cigar = (char *)malloc((strlen(cigar)+1)* sizeof(char));
        strcpy(read[i].cigar,cigar);

        read[i].length = strlen(seq);
        read[i].seq = (char *)malloc((read[i].length+1)*sizeof(char));
        read[i].qual = (char *)malloc((strlen(qual)+1)* sizeof(char));
        strcpy(read[i].seq,seq);
        strcpy(read[i].qual,qual);

        read[i].NM = 0;
        read[i].score = 0;

        read[i].site--;

        i++;
    }
    free(f_line);
    free(name);
    free(cigar);
    free(seq);
    free(qual);
    return i;
}
void free_read_sam(struct read_inf_t *read,int read_num)
{
    int i;
    for(i = 0;i<read_num;i++)
    {
        if (read->name!=NULL){free(read->name);read->name = NULL;}
        if (read->cigar!=NULL){free(read->cigar);read->cigar = NULL;}
        if (read->seq!=NULL) {free(read->seq);read->seq = NULL;}
        if (read->qual!=NULL) {free(read->qual);read->qual = NULL;}
    }
}
void output_read(struct read_inf_t *read,int read_num,FILE *out)//,FILE *un)
{
    int i;
    for(i = 0;i<read_num;i++)
    {
        if(read[i].chr!=-1)
        fprintf(out,"%s\t%d\t%s\t%lu\t0\t%s\t*\t0\t0\t%s\t%s\t%s\n",
            read[i].name,read[i].flag,opt->chr->list[read[i].chr].name,read[i].site+1,read[i].cigar,read[i].seq,read[i].qual,read[i].TAG);
        else
        fprintf(out,"%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\t%s\n",
            read[i].name,read[i].flag,read[i].seq,read[i].qual,read[i].TAG);
    }
}
struct exon_t
{
    unsigned int start;
    unsigned int end;
    char name[100];
};

int load_gtf(struct m_opt *opt)
{
    char *f_line = (char *)malloc(MAX_STRING_LENGTH*sizeof(char));
    char chr[100];
    char name[100];
    unsigned int start;
    unsigned int end;
    int chr_order = 0;
    unsigned int chr_start = 0;

    struct exon_t exon;

    struct splice_list *SP = (struct splice_list *)malloc(opt->SP_MAX*sizeof(struct splice_list));
	int SP_num = 0;
    struct splice_list sp;

    FILE *file;
    file = fopen(opt->gtf_file,"r");

    while(fgets(f_line,MAX_STRING_LENGTH,file)!=NULL)
    {
        sscanf(f_line,"%s\t%s%*s\t%u\t%u",chr,name,&start,&end);
        chr_order = find_chr_name(chr);
        if(chr_order<0||chr_order>=opt->chr->total) return 1;

        chr_start = opt->chr->list[chr_order].start_site;
        start += chr_start;
        end += chr_start;

        if(strcmp(name,exon.name)==0)
        {
            if(exon.end < start)
            {
                sp.start = exon.end+1;
                sp.end = start-1;
                sp.num = 1000;
                Update_SP_s(SP,&SP_num,opt->SP_MAX,&sp,0);
            }
        }
        else
            strcpy(exon.name,name);

        exon.start = start;
        exon.end = end;
    }
    free(SP);
    free(f_line);
    return 0;
}
void process_cigar(char *cigarBuf,struct cigar_t *cigar,int *cigar_total)//,unsigned int *end_pos)
{
    int i,j;
    char temp[20];
    {
        j = 0;
        (*cigar_total) = 0;
        for(i = 0;i<strlen(cigarBuf);i++)
        {
            if((cigarBuf[i]>='0')&&(cigarBuf[i]<='9'))
            {
                temp[j] = cigarBuf[i];
                j++;
            }
            else
            {
                temp[j] = '\0';
                cigar[(*cigar_total)].c = cigarBuf[i];
                cigar[(*cigar_total)].l = atoi(temp);

                //if((cigar[(*cigar_total)].c=='M')||(cigar[(*cigar_total)].c=='X')||(cigar[(*cigar_total)].c=='D')||(cigar[(*cigar_total)].c=='N')) (*end_pos)+= cigar[(*cigar_total)].l;
                (*cigar_total)++;
                j = 0;
            }
        }
    }
}
int re_check_exon(struct read_inf_t *read0,char *read,int length,struct cigar_t *cigar,int *cigar_num,uint64_t start_pos,int chr_order,struct Splice_DP_t *DP,char *cigarS)
{
    int break_flag = 1;
    int EXON_LENGTH = 0;//int EXON_LENGTH = 50;

    int i,j;
    unsigned int ref_pos = start_pos;
    int read_pos = 0;

    struct exon_seed_t *seed = (struct exon_seed_t *)malloc(EXON_SEED*sizeof(struct exon_seed_t));
    int seed_num = 0;
    int n = 0;
    int score = 0;
    int l = 0;

    struct cigar_t *t_cigar = (struct cigar_t *)malloc(MAX_CIGAR_BUF*sizeof(struct cigar_t));
    int t_num = 0;

    int seed_order[1000];
    int seed_order_n = 0;

    int max_s,max_o;

    seed[0].ref_start = start_pos;
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

    max_o = seed_num-1;
    seed_order_n = seed_num;
    for(i = 0;i<seed_num;i++)
        seed_order[i]=i;

    //check_splice;
    int signal_flag = 0;//0 GT-AG 1 CT-AC 2 unusual
    int num_1 = 0,num_2 = 0;
    int order[300][3];
    int el;
    int k;

    seed[seed_order[0]].c_start = 0;
    seed[seed_order[seed_order_n-1]].c_end = seed[seed_order[seed_order_n-1]].cigar_num-1;
    for(i = 1;i<seed_order_n;i++)//head
    {
        el = 0;
        for(j = 0;j<seed[seed_order[i]].cigar_num;j++)
        {
            if((seed[seed_order[i]].cigar[j].c=='M')&&(seed[seed_order[i]].cigar[j].l>=8)&&(el+seed[seed_order[i]].cigar[j].l>=10)) break;

            if((seed[seed_order[i]].cigar[j].c=='M')||(seed[seed_order[i]].cigar[j].c=='X')||(seed[seed_order[i]].cigar[j].c=='I')||(seed[seed_order[i]].cigar[j].c=='S'))el+=seed[seed_order[i]].cigar[j].l;
        }
        if(j<seed[seed_order[i]].cigar_num)
        {
            k = j;
        j = 0;
        while(j<=k)
        {
            if(j==k)
            {
                seed[seed_order[i]].ref_start+= 4;
                seed[seed_order[i]].read_start+= 4;
                break;
            }
            if(seed[seed_order[i]].cigar[j].c=='M') {seed[seed_order[i]].ref_start+=seed[seed_order[i]].cigar[j].l;seed[seed_order[i]].read_start+=seed[seed_order[i]].cigar[j].l;}
            else if(seed[seed_order[i]].cigar[j].c=='X'){seed[seed_order[i]].ref_start+=seed[seed_order[i]].cigar[j].l;seed[seed_order[i]].read_start+=seed[seed_order[i]].cigar[j].l;}
            else if(seed[seed_order[i]].cigar[j].c=='I'){seed[seed_order[i]].read_start+=seed[seed_order[i]].cigar[j].l;}
            else if(seed[seed_order[i]].cigar[j].c=='S'){seed[seed_order[i]].read_start+=seed[seed_order[i]].cigar[j].l;}
            else if(seed[seed_order[i]].cigar[j].c=='D'){seed[seed_order[i]].ref_start+=seed[seed_order[i]].cigar[j].l;}
            j++;
        }
        seed[seed_order[i]].c_start = j;
        }
        else
        {
            max_s = 0;
            max_o = 0;
            for(j = 0;j<seed[seed_order[i]].cigar_num;j++)
            {
                if(seed[seed_order[i]].cigar[j].c=='M')
                {
                    if(seed[seed_order[i]].cigar[j].l>max_s)
                    {
                        max_s = seed[seed_order[i]].cigar[j].l;
                        max_o = j;
                    }
                }
            }
            for(j = 0;j<max_o;j++)
            {
            if((seed[seed_order[i]].cigar[j].c=='M')||(seed[seed_order[i]].cigar[j].c=='X')){seed[seed_order[i]].ref_start+=seed[seed_order[i]].cigar[j].l;seed[seed_order[i]].read_start+=seed[seed_order[i]].cigar[j].l;}
            else if(seed[seed_order[i]].cigar[j].c=='I'){seed[seed_order[i]].read_start+=seed[seed_order[i]].cigar[j].l;}
            else if(seed[seed_order[i]].cigar[j].c=='S'){seed[seed_order[i]].read_start+=seed[seed_order[i]].cigar[j].l;}
            else if(seed[seed_order[i]].cigar[j].c=='D'){seed[seed_order[i]].ref_start+=seed[seed_order[i]].cigar[j].l;}
            }
            seed[seed_order[i]].c_start = max_o;
            seed[seed_order[i]].ref_start+= 4;
            seed[seed_order[i]].read_start+= 4;
        }
    }
    for(i = 0;i<seed_order_n-1;i++)//tail
    {
        el = 0;
        for(j = seed[seed_order[i]].cigar_num-1;j>=seed[seed_order[i]].c_start;j--)
        {
            if((seed[seed_order[i]].cigar[j].c=='M')&&(seed[seed_order[i]].cigar[j].l>=8)&&(el+seed[seed_order[i]].cigar[j].l>10)) break;

            if((seed[seed_order[i]].cigar[j].c=='M')||(seed[seed_order[i]].cigar[j].c=='X')||(seed[seed_order[i]].cigar[j].c=='I')||(seed[seed_order[i]].cigar[j].c=='S'))el+=seed[seed_order[i]].cigar[j].l;
        }
        if(j>=seed[seed_order[i]].c_start)
        {
        k = j;
        j = seed[seed_order[i]].cigar_num-1;
        while(j>=k)
        {
            if(j==k)
            {
                seed[seed_order[i]].cigar[j].l-=4;
                seed[seed_order[i]].ref_end-= 4;
                seed[seed_order[i]].read_end-= 4;
                break;
            }
            if(seed[seed_order[i]].cigar[j].c=='M') {seed[seed_order[i]].ref_end-=seed[seed_order[i]].cigar[j].l;seed[seed_order[i]].read_end-=seed[seed_order[i]].cigar[j].l;}
            else if(seed[seed_order[i]].cigar[j].c=='X'){seed[seed_order[i]].ref_end-=seed[seed_order[i]].cigar[j].l;seed[seed_order[i]].read_end-=seed[seed_order[i]].cigar[j].l;}
            else if(seed[seed_order[i]].cigar[j].c=='I'){seed[seed_order[i]].read_end-=seed[seed_order[i]].cigar[j].l;}
            else if(seed[seed_order[i]].cigar[j].c=='S'){seed[seed_order[i]].read_end-=seed[seed_order[i]].cigar[j].l;}
            else if(seed[seed_order[i]].cigar[j].c=='D'){seed[seed_order[i]].ref_end-=seed[seed_order[i]].cigar[j].l;}
            j--;
        }
        seed[seed_order[i]].c_end = k;
        }
        else
        {
            max_s = 0;
            max_o = 0;
            for(j = seed[seed_order[i]].cigar_num-1;j>=seed[seed_order[i]].c_start;j--)
            {
                if(seed[seed_order[i]].cigar[j].c=='M')
                {
                    if(seed[seed_order[i]].cigar[j].l>max_s)
                    {
                        max_s = seed[seed_order[i]].cigar[j].l;
                        max_o = j;
                    }
                }
            }
            for(j = seed[seed_order[i]].cigar_num-1;j>max_o;j--)
            {
            if((seed[seed_order[i]].cigar[j].c=='M')||(seed[seed_order[i]].cigar[j].c=='X')){seed[seed_order[i]].ref_end-=seed[seed_order[i]].cigar[j].l;seed[seed_order[i]].read_end-=seed[seed_order[i]].cigar[j].l;}
            else if(seed[seed_order[i]].cigar[j].c=='I'){seed[seed_order[i]].read_end-=seed[seed_order[i]].cigar[j].l;}
            else if(seed[seed_order[i]].cigar[j].c=='S'){seed[seed_order[i]].read_end-=seed[seed_order[i]].cigar[j].l;}
            else if(seed[seed_order[i]].cigar[j].c=='D'){seed[seed_order[i]].ref_end-=seed[seed_order[i]].cigar[j].l;}
            }

            seed[seed_order[i]].c_end = j;
            if (j==seed[seed_order[i]].c_start)
            {
            if(max_s>=8)
            {
                seed[seed_order[i]].cigar[j].l-=4;
                seed[seed_order[i]].ref_end-= 4;
                seed[seed_order[i]].read_end-= 4;
            }
            else if(max_s>=4)
            {
                seed[seed_order[i]].cigar[j].l-=max_s-4;
                seed[seed_order[i]].ref_end-= max_s-4;
                seed[seed_order[i]].read_end-= max_s-4;
            }
            else
            {
                seed[seed_order[i]].cigar[j].l-=0;
            }
            }
            else
            {
                max_o=min(max_s,4);
                seed[seed_order[i]].cigar[j].l-=max_o;
                seed[seed_order[i]].ref_end-= max_o;
                seed[seed_order[i]].read_end-= max_o;
            }
        }
    }
    for(i = 1;i<seed_order_n;i++)
    {
        seed[seed_order[i]].cigar[seed[seed_order[i]].c_start].l-=4;
        if(seed[seed_order[i]].cigar[seed[seed_order[i]].c_start].l<0)
        {
            seed[seed_order[i]].ref_start+= seed[seed_order[i]].cigar[seed[seed_order[i]].c_start].l;
            seed[seed_order[i]].read_start+= seed[seed_order[i]].cigar[seed[seed_order[i]].c_start].l;
            seed[seed_order[i]].cigar[seed[seed_order[i]].c_start].l=0;
        }
    }

    for(i = 1;i<seed_order_n;i++)
    {
        order[i-1][0] = -1;
        order[i-1][1] = -1;
        order[i-1][2] = -1;

        max_s = -1;

        seed[seed_num].ref_start = seed[seed_order[i-1]].ref_end+1;
        seed[seed_num].ref_end = seed[seed_order[i]].ref_start-1;
        seed[seed_num].read_start = seed[seed_order[i-1]].read_end+1;
        seed[seed_num].read_end = seed[seed_order[i]].read_start-1;

        seed[seed_num].c_start = 0;
        seed[seed_num].c_end = 0;
        seed[seed_num].score = 0;

        l = seed[seed_num].read_end-seed[seed_num].read_start+1;
        t_num = 0;
        splice_pos(chr_order,seed[seed_num].ref_start,seed[seed_num].ref_end,read+seed[seed_num].read_start,l,t_cigar,&t_num,&(seed[seed_num].score),MAX_CIGAR_BUF-1,0,DP);
        if(seed[seed_num].score!=-1000)
        {
            seed[seed_num].score+=16;
            seed[seed_num].cigar = (struct cigar_t *)malloc(t_num*sizeof(struct cigar_t));
            memcpy(seed[seed_num].cigar,t_cigar,t_num*sizeof(struct cigar_t));
            seed[seed_num].cigar_num = t_num;
            order[i-1][0] = seed_num;
            max_s = seed_num;
            seed_num++;
            if(seed_num>=EXON_SEED)
            {
            for(i = 0;i<seed_num;i++)
            {
                if(seed[i].cigar!=NULL)free(seed[i].cigar);
            }
            free(seed);
            free(t_cigar);
            cigarS[0]='\0';
            return 1;
            }
        }

        seed[seed_num].ref_start = seed[seed_order[i-1]].ref_end+1;
        seed[seed_num].ref_end = seed[seed_order[i]].ref_start-1;
        seed[seed_num].read_start = seed[seed_order[i-1]].read_end+1;
        seed[seed_num].read_end = seed[seed_order[i]].read_start-1;

        seed[seed_num].c_start = 0;
        seed[seed_num].c_end = 0;
        seed[seed_num].score = 0;

        t_num = 0;
        splice_pos(chr_order,seed[seed_num].ref_start,seed[seed_num].ref_end,read+seed[seed_num].read_start,l,t_cigar,&t_num,&(seed[seed_num].score),MAX_CIGAR_BUF-1,1,DP);
        if(seed[seed_num].score!=-1000)
        {
            seed[seed_num].score+=16;
            seed[seed_num].cigar = (struct cigar_t *)malloc(t_num*sizeof(struct cigar_t));
            memcpy(seed[seed_num].cigar,t_cigar,t_num*sizeof(struct cigar_t));
            seed[seed_num].cigar_num = t_num;
            order[i-1][1] = seed_num;
            if(max_s==-1) max_s = seed_num;
            else if(seed[seed_num].score>seed[seed_num-1].score) max_s = seed_num;

            seed_num++;
            if(seed_num>=EXON_SEED)
            {
            for(i = 0;i<seed_num;i++)
            {
                if(seed[i].cigar!=NULL)free(seed[i].cigar);
            }
            free(seed);
            free(t_cigar);
            cigarS[0]='\0';
            return 1;
            }
        }

        seed[seed_num].ref_start = seed[seed_order[i-1]].ref_end+1;
        seed[seed_num].ref_end = seed[seed_order[i]].ref_start-1;
        seed[seed_num].read_start = seed[seed_order[i-1]].read_end+1;
        seed[seed_num].read_end = seed[seed_order[i]].read_start-1;

        seed[seed_num].c_start = 0;
        seed[seed_num].c_end = 0;
        seed[seed_num].score = 0;

        t_num = 0;
        if(seed[seed_num].ref_end-seed[seed_num].ref_start+1>max(l*(1.1),l+20))
        seed[seed_num].score = Splice_DP(opt->chr->list[chr_order].seq,seed[seed_num].ref_start,seed[seed_num].ref_end,read+seed[seed_num].read_start,l,t_cigar,&t_num,DP);
        else
        seed[seed_num].score = Area_DP(opt->chr->list[chr_order].seq,seed[seed_num].ref_start,seed[seed_num].ref_end,read+seed[seed_num].read_start,l,t_cigar,&t_num,DP);

        if(seed[seed_num].score!=-10000)
        {
            seed[seed_num].score = cigar_score(t_cigar,t_num)+opt->splice;
            seed[seed_num].cigar = (struct cigar_t *)malloc(t_num*sizeof(struct cigar_t));
            memcpy(seed[seed_num].cigar,t_cigar,t_num*sizeof(struct cigar_t));
            seed[seed_num].cigar_num = t_num;

            order[i-1][2] = seed_num;
            seed_num++;
        }
        if(seed_num>=EXON_SEED)
        {
            for(i = 0;i<seed_num;i++)
            {
                if(seed[i].cigar!=NULL)free(seed[i].cigar);
            }
            free(seed);
            free(t_cigar);
            cigarS[0]='\0';
            return 1;
        }

        seed[seed_num].ref_start = seed[seed_order[i-1]].ref_end+1;
        seed[seed_num].ref_end = seed[seed_order[i]].ref_start-1;
        seed[seed_num].read_start = seed[seed_order[i-1]].read_end+1;
        seed[seed_num].read_end = seed[seed_order[i]].read_start-1;

        seed[seed_num].c_start = 0;
        seed[seed_num].c_end = 0;
        seed[seed_num].score = 0;

        t_num = 0;
        Known_splice_pos(chr_order,seed[seed_num].ref_start,seed[seed_num].ref_end,read+seed[seed_num].read_start,l,t_cigar,&t_num,&(seed[seed_num].score),MAX_CIGAR_BUF-1,DP);
        if(seed[seed_num].score!=-1000)
        {
            seed[seed_num].score+=16;
            if(seed[seed_num].score<seed[seed_num-1].score) seed_num--;
            else
            {
                seed[seed_num].cigar = (struct cigar_t *)malloc(t_num*sizeof(struct cigar_t));
                memcpy(seed[seed_num].cigar,t_cigar,t_num*sizeof(struct cigar_t));
                seed[seed_num].cigar_num = t_num;
            }

            order[i-1][2] = seed_num;
            if(max_s==-1) max_s = seed_num;
            else if(seed[seed_num].score>seed[max_s].score) max_s = seed_num;
            seed_num++;
        }
        if(seed_num>=EXON_SEED)
        {
            for(i = 0;i<seed_num;i++)
            {
                if(seed[i].cigar!=NULL)free(seed[i].cigar);
            }
            free(seed);
            free(t_cigar);
            cigarS[0]='\0';
            return 1;
        }

        if(order[i-1][0]==max_s) num_1++;
        else if(order[i-1][1]==max_s) num_2++;
    }
    //write cigar
    char temp[20];
    cigarS[0] = '\0';
    int M = 0;

    //int lengthS = 0;

    unsigned int start = 0;

    if(num_1>num_2)signal_flag = 0;//GT-AG
    else if(num_1==num_2) signal_flag = 2;
    else signal_flag = 1;
    {
    int f_seed_order[1000];
    int f_seed_order_n = 0;

    for(i = 0;i<seed_order_n;i++)
    {
        if(seed[seed_order[i]].cigar[seed[seed_order[i]].c_start].l==0) seed[seed_order[i]].c_start++;
        if(seed[seed_order[i]].cigar[seed[seed_order[i]].c_end].l==0) seed[seed_order[i]].c_end--;
        f_seed_order[f_seed_order_n] = seed_order[i];
        f_seed_order_n++;

        if(i!=seed_order_n-1)
        {
            if(signal_flag==2)
            {
                max_s = order[i][0];
                if(max_s==-1) max_s = order[i][1];
                else if((max_s!=-1)&&(order[i][1]!=-1)&&(seed[order[i][1]].score>seed[max_s].score)) max_s = order[i][1];
            }
            else max_s = order[i][signal_flag];

            if(max_s==-1) max_s = order[i][2];
            else if((max_s!=-1)&&(order[i][2]!=-1)&&(seed[order[i][2]].score>seed[max_s].score)) max_s = order[i][2];

            if(max_s==-1)
            {
            for(i = 0;i<seed_num;i++)
            {
                if(seed[i].cigar!=NULL)free(seed[i].cigar);
            }
            free(seed);
            free(t_cigar);
            cigarS[0]='\0';
            return 1;
            }

            seed[max_s].c_start = 0;
            seed[max_s].c_end = seed[max_s].cigar_num-1;

            f_seed_order[f_seed_order_n] = max_s;
            f_seed_order_n++;
        }
        if(seed[seed_order[i]].read_end<seed[seed_order[i]].read_start)//0M
        {
            if(i==0)
            {
                f_seed_order_n = 1;
                f_seed_order[0] = max_s;
                seed[f_seed_order[f_seed_order_n-1]].c_start++;
            }
            else if(i==seed_order_n-1)
            {
                f_seed_order_n--;
                if (seed[f_seed_order[f_seed_order_n-1]].cigar[seed[f_seed_order[f_seed_order_n-1]].c_end].c=='N')
                    seed[f_seed_order[f_seed_order_n-1]].c_end--;
            }
            else
            {
                if((f_seed_order_n>=3)&&(seed[f_seed_order[f_seed_order_n-1]].cigar[seed[f_seed_order[f_seed_order_n-1]].c_start].c=='N')&&(seed[f_seed_order[f_seed_order_n-3]].cigar[seed[f_seed_order[f_seed_order_n-3]].c_end].c=='N'))
        {
        order[seed_order_n][0] = -1;
        order[seed_order_n][1] = -1;
        order[seed_order_n][2] = -1;

        max_s = -1;

        seed[seed_num].ref_start = seed[seed_order[i-1]].ref_end+1;
        seed[seed_num].ref_end = seed[seed_order[i+1]].ref_start-1;
        seed[seed_num].read_start = seed[seed_order[i-1]].read_end+1;
        seed[seed_num].read_end = seed[seed_order[i+1]].read_start-1;

        seed[seed_num].c_start = 0;
        seed[seed_num].c_end = 0;
        seed[seed_num].score = 0;

        l = seed[seed_num].read_end-seed[seed_num].read_start+1;
        if(l<0)
            printf("1\n");
        t_num = 0;
        splice_pos(chr_order,seed[seed_num].ref_start,seed[seed_num].ref_end,read+seed[seed_num].read_start,l,t_cigar,&t_num,&(seed[seed_num].score),MAX_CIGAR_BUF-1,0,DP);
        if(seed[seed_num].score!=-1000)
        {
            seed[seed_num].score+=16;
            seed[seed_num].cigar = (struct cigar_t *)malloc(t_num*sizeof(struct cigar_t));
            memcpy(seed[seed_num].cigar,t_cigar,t_num*sizeof(struct cigar_t));
            seed[seed_num].cigar_num = t_num;
            order[seed_order_n][0] = seed_num;
            max_s = seed_num;
            seed_num++;
            if(seed_num>=EXON_SEED)
            {
            for(i = 0;i<seed_num;i++)
            {
                if(seed[i].cigar!=NULL)free(seed[i].cigar);
            }
            free(seed);
            free(t_cigar);
            cigarS[0]='\0';
            return 1;
            }
        }

        seed[seed_num].ref_start = seed[seed_order[i-1]].ref_end+1;
        seed[seed_num].ref_end = seed[seed_order[i+1]].ref_start-1;
        seed[seed_num].read_start = seed[seed_order[i-1]].read_end+1;
        seed[seed_num].read_end = seed[seed_order[i+1]].read_start-1;

        seed[seed_num].c_start = 0;
        seed[seed_num].c_end = 0;
        seed[seed_num].score = 0;

        t_num = 0;
        splice_pos(chr_order,seed[seed_num].ref_start,seed[seed_num].ref_end,read+seed[seed_num].read_start,l,t_cigar,&t_num,&(seed[seed_num].score),MAX_CIGAR_BUF-1,1,DP);
        if(seed[seed_num].score!=-1000)
        {
            seed[seed_num].score+=16;
            seed[seed_num].cigar = (struct cigar_t *)malloc(t_num*sizeof(struct cigar_t));
            memcpy(seed[seed_num].cigar,t_cigar,t_num*sizeof(struct cigar_t));
            seed[seed_num].cigar_num = t_num;
            order[seed_order_n][1] = seed_num;
            if(max_s==-1) max_s = seed_num;
            else if(seed[seed_num].score>seed[seed_num-1].score) max_s = seed_num;

            seed_num++;
            if(seed_num>=EXON_SEED)
            {
            for(i = 0;i<seed_num;i++)
            {
                if(seed[i].cigar!=NULL)free(seed[i].cigar);
            }
            free(seed);
            free(t_cigar);
            cigarS[0]='\0';
            return 1;
            }
        }

        seed[seed_num].ref_start = seed[seed_order[i-1]].ref_end+1;
        seed[seed_num].ref_end = seed[seed_order[i+1]].ref_start-1;
        seed[seed_num].read_start = seed[seed_order[i-1]].read_end+1;
        seed[seed_num].read_end = seed[seed_order[i+1]].read_start-1;

        seed[seed_num].c_start = 0;
        seed[seed_num].c_end = 0;
        seed[seed_num].score = 0;

        t_num = 0;
        if(seed[seed_num].ref_end-seed[seed_num].ref_start+1>2*max(l*(1.1),l+20))
        seed[seed_num].score = Splice_DP(opt->chr->list[chr_order].seq,seed[seed_num].ref_start,seed[seed_num].ref_end,read+seed[seed_num].read_start,l,t_cigar,&t_num,DP);
        else
        seed[seed_num].score = Area_DP(opt->chr->list[chr_order].seq,seed[seed_num].ref_start,seed[seed_num].ref_end,read+seed[seed_num].read_start,l,t_cigar,&t_num,DP);

        if(seed[seed_num].score!=-10000)
        {
            seed[seed_num].score = cigar_score(t_cigar,t_num)+opt->splice;
            seed[seed_num].cigar = (struct cigar_t *)malloc(t_num*sizeof(struct cigar_t));
            memcpy(seed[seed_num].cigar,t_cigar,t_num*sizeof(struct cigar_t));
            seed[seed_num].cigar_num = t_num;

            order[seed_order_n][2] = seed_num;
            seed_num++;
        }
        if(seed_num>=EXON_SEED)
        {
            for(i = 0;i<seed_num;i++)
            {
                if(seed[i].cigar!=NULL)free(seed[i].cigar);
            }
            free(seed);
            free(t_cigar);
            cigarS[0]='\0';
            return 1;
        }

        seed[seed_num].ref_start = seed[seed_order[i-1]].ref_end+1;
        seed[seed_num].ref_end = seed[seed_order[i+1]].ref_start-1;
        seed[seed_num].read_start = seed[seed_order[i-1]].read_end+1;
        seed[seed_num].read_end = seed[seed_order[i+1]].read_start-1;

        seed[seed_num].c_start = 0;
        seed[seed_num].c_end = 0;
        seed[seed_num].score = 0;

        t_num = 0;
        Known_splice_pos(chr_order,seed[seed_num].ref_start,seed[seed_num].ref_end,read+seed[seed_num].read_start,l,t_cigar,&t_num,&(seed[seed_num].score),MAX_CIGAR_BUF-1,DP);
        if(seed[seed_num].score!=-1000)
        {
            seed[seed_num].score+=16;
            if(seed[seed_num].score<seed[seed_num-1].score) seed_num--;
            else
            {
                seed[seed_num].cigar = (struct cigar_t *)malloc(t_num*sizeof(struct cigar_t));
                memcpy(seed[seed_num].cigar,t_cigar,t_num*sizeof(struct cigar_t));
                seed[seed_num].cigar_num = t_num;
            }
            order[seed_order_n][2] = seed_num;
            if(max_s==-1) max_s = seed_num;
            else if(seed[seed_num].score>seed[max_s].score) max_s = seed_num;
            seed_num++;
        }
        if(seed_num>=EXON_SEED)
        {
            for(i = 0;i<seed_num;i++)
            {
                if(seed[i].cigar!=NULL)free(seed[i].cigar);
            }
            free(seed);
            free(t_cigar);
            cigarS[0]='\0';
            return 1;
        }

        f_seed_order_n-=3;

            if(signal_flag==2)
            {
                max_s = order[seed_order_n][0];
                if(max_s==-1) max_s = order[seed_order_n][1];
                else if((max_s!=-1)&&(order[seed_order_n][1]!=-1)&&(seed[order[seed_order_n][1]].score>seed[max_s].score)) max_s = order[seed_order_n][1];
            }
            else max_s = order[seed_order_n][signal_flag];

            if(max_s==-1) max_s = order[seed_order_n][2];
            else if((max_s!=-1)&&(order[seed_order_n][2]!=-1)&&(seed[order[seed_order_n][2]].score>seed[max_s].score)) max_s = order[seed_order_n][2];

            if(max_s==-1)
            {
            for(i = 0;i<seed_num;i++)
            {
                if(seed[i].cigar!=NULL)free(seed[i].cigar);
            }
            free(seed);
            free(t_cigar);
            cigarS[0]='\0';
            return 1;
            }

            seed[max_s].c_start = 0;
            seed[max_s].c_end = seed[max_s].cigar_num-1;

            f_seed_order[f_seed_order_n] = max_s;
            f_seed_order_n++;
    }
            else
                {
                    f_seed_order[f_seed_order_n-2] = f_seed_order[f_seed_order_n-1];
                    f_seed_order_n--;
                }
            }
        }
    }
    seed_order_n = f_seed_order_n;
    memcpy(seed_order,f_seed_order,seed_order_n*sizeof(int));

    for(i = 0;i<seed_order_n;i++)
    {
        //check_cigar2(seed[seed_order[i]].cigar+seed[seed_order[i]].c_start,seed[seed_order[i]].c_end-seed[seed_order[i]].c_start+1,seed[seed_order[i]].read_end-seed[seed_order[i]].read_start+1);
        start = seed[seed_order[i]].ref_start;
        for(j = seed[seed_order[i]].c_start;j<=seed[seed_order[i]].c_end;j++)
        {
                if((seed[seed_order[i]].cigar[j].c=='M')||(seed[seed_order[i]].cigar[j].c=='X')) {M+=seed[seed_order[i]].cigar[j].l;start+=seed[seed_order[i]].cigar[j].l;}
                else
                {
                    if(M!=0)
                    {
                        sprintf(temp,"%dM",M);
                        strcat(cigarS,temp);
                        M=0;
                    }
                    if((seed[seed_order[i]].cigar[j].l<=8)&&(seed[seed_order[i]].cigar[j].c=='N')) seed[seed_order[i]].cigar[j].c = 'D';
                    if((seed[seed_order[i]].cigar[j].l>8)&&(seed[seed_order[i]].cigar[j].c=='D')) seed[seed_order[i]].cigar[j].c = 'N';
                    sprintf(temp,"%d%c",seed[seed_order[i]].cigar[j].l,seed[seed_order[i]].cigar[j].c);
                    strcat(cigarS,temp);
                    if((seed[seed_order[i]].cigar[j].c == 'N')||(seed[seed_order[i]].cigar[j].c == 'D')) start+=seed[seed_order[i]].cigar[j].l;
                }

                if(seed[seed_order[i]].cigar[j].c=='M') read0->score+=seed[seed_order[i]].cigar[j].l*opt->match;
                else if(seed[seed_order[i]].cigar[j].c=='X') {read0->score-=seed[seed_order[i]].cigar[j].l*opt->miss; read0->NM+=seed[seed_order[i]].cigar[j].l;}
                else if(seed[seed_order[i]].cigar[j].c=='I') {read0->score-= 4+(seed[seed_order[i]].cigar[j].l-1)*opt->gap; read0->NM+=seed[seed_order[i]].cigar[j].l;}
                else if(seed[seed_order[i]].cigar[j].c=='D') read0->score-= 4+(seed[seed_order[i]].cigar[j].l-1)*opt->gap;
                else if(seed[seed_order[i]].cigar[j].c=='N') read0->score-= opt->splice;
        }
    }
    }

    if(M!=0)
    {
        sprintf(temp,"%dM",M);
        strcat(cigarS,temp);
    }

    pthread_mutex_lock(&UpdateLock);
    qsort(opt->SP,opt->SP_num,sizeof(struct splice_list),SP_cmp);
    pthread_mutex_unlock(&UpdateLock);

    for(i = 0;i<seed_num;i++)
    {
        if(seed[i].cigar!=NULL)free(seed[i].cigar);
    }
    free(seed);
    free(t_cigar);
    return 0;
}
void change_cigar(struct read_inf_t *read0,struct cigar_t *cigar,int cigar_num,char *cigarS)
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
            if((cigar[i].l<=8)&&(cigar[i].c=='N')) cigar[i].c = 'D';
            if((cigar[i].l>8)&&(cigar[i].c=='D')) cigar[i].c = 'N';

            sprintf(temp,"%d%c",cigar[i].l,cigar[i].c);
            strcat(cigarS,temp);
        }

        if(cigar[i].c=='M') read0->score+=cigar[i].l*opt->match;
        else if(cigar[i].c=='X') {read0->score-=cigar[i].l*opt->miss; read0->NM+=cigar[i].l;}
        else if(cigar[i].c=='I') {read0->score-= 4+(cigar[i].l-1)*opt->gap; read0->NM+=cigar[i].l;}
        else if(cigar[i].c=='D') read0->score-= 4+(cigar[i].l-1)*opt->gap;
        else if(cigar[i].c=='N') read0->score-= opt->splice;
    }
    if(M!=0)
    {
        sprintf(temp,"%dM",M);
        strcat(cigarS,temp);
    }
}
void *recheck_core(void* arg)
{
    struct job_seed *job = (struct job_seed *)arg;

	struct read_inf_t *read = (struct read_inf_t *)malloc(opt->thread_block*sizeof(struct read_inf_t));
	int read_num = 0;

	struct cigar_t *cigar = (struct cigar_t *)calloc(MAX_CIGAR_BUF,sizeof(struct cigar_t));
	int cigar_num;

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

	char *cigarS = (char *)malloc((MAX_STRING_LENGTH*2)*sizeof(char));
	int lengthS;

	while (1)
	{
        //read_file;
        pthread_mutex_lock(&ReadLock);
        read_num = read_file_sam(read,job->in_file);
		pthread_mutex_unlock(&ReadLock);

		if(read_num<=0) break;
        for (i = 0;i<read_num;i++)
        {
            if(read[i].cigar[0]=='*')
            {
                sprintf(read[i].TAG,"NM:i:0\tAS:i:0");
                continue;
            }
            cigar_num = 0;
            process_cigar(read[i].cigar,cigar,&cigar_num);

            cigarS[0] = '\0';
            if(re_check_exon(read+i,read[i].seq,read[i].length,cigar,&cigar_num,read[i].site,read[i].chr,DP,cigarS)==1)
            {
                fprintf(stdout,"%s check exon fail\n",read[i].name);
                change_cigar(read+i,cigar,cigar_num,cigarS);
            }

            lengthS = strlen(cigarS);
            free(read[i].cigar);
            read[i].cigar = (char *)calloc(lengthS+1,sizeof(char));
            strcpy(read[i].cigar,cigarS);
            read[i].cigar[lengthS] = '\0';

            sprintf(read[i].TAG,"NM:i:%d\tAS:i:%d",read[i].NM,read[i].score);
        }

        pthread_mutex_lock(&OutputLock);
        output_read(read,read_num,job->out_file);
        pthread_mutex_unlock(&OutputLock);

        free_read_sam(read,read_num);
	}

	free(read);
	free(cigar);
	free(cigarS);

    for (i = 0; i <=MAX_DP_AREA; i++)
	{
		free(DP->r[i]); free(DP->t[i]);free(DP->r1[i]); free(DP->t1[i]); free(DP->d[i]); free(DP->a[i]);
	}
	free(DP->r); free(DP->t);free(DP->r1); free(DP->t1); free(DP->d); free(DP->a);
	free(DP->cigarD);free(DP-> max_d);free(DP-> max_a);
	free(DP);

    return (void*)(0);
}
void out_put_SAM_head(struct m_opt *opt,FILE *file)
{
    fprintf(file,"@HD\tVN:1.0\tSO:unsorted\n");
    fprintf(file,"@PG\tID:DEEPL\tVN:2020-07-01\n");
    int i = 0;
    for(i = 0;i<opt->chr->total;i++)
        fprintf(file,"@SQ\tSN:%s\tLN:%u\n",opt->chr->list[i].name,opt->chr->list[i].length);
}
int re_check(struct m_opt *opt)
{
    struct job_seed *job = (struct job_seed *)calloc(1,sizeof(struct job_seed));
    int i;

    char temp_name[MAX_NAME_LENGTH+50];
    sprintf(temp_name,"%s.temp",opt->Output_path);
    job->in_file = fopen(temp_name,"r");
    if (job->in_file == NULL) return 1;

    job->out_file = fopen(opt->Output_path,"w");
    if (job->out_file == NULL){return 1;}

    out_put_SAM_head(opt,job->out_file);

	pthread_t *pthreads = malloc(sizeof(pthread_t) * opt->thread_num);

    for (i = 0; i < opt->thread_num; i++) pthread_create(&pthreads[i], NULL, recheck_core, job);
	for (i = 0; i < opt->thread_num; i++) pthread_join(pthreads[i], NULL);

    fclose(job->in_file);
    remove(temp_name);
    fclose(job->out_file);

    free(pthreads);
    free(job);

    return 0;
}
