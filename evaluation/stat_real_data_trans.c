#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAX_NAME_LENGTH 400
#define READ_MAX_LENGTH 400
#define MAX_STRING_LENGTH 40000

#define MAX_EXON 5000000
#define MAX_TRANS 1000000

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

struct exon_t
{
    char chr[40];
    unsigned int start;
    unsigned int end;
    int flag;
    int exact;
    int part;
};
struct trans_t
{
    char name[100];
    char chr[40];
    unsigned int start;
    unsigned int end;
    int flag;

    struct exon_t *exon;
    int exon_num;
};
int find_exon(struct exon_t *SP,int num,char *chr,unsigned int start)
{
    int left = 0, right = num-1, middle;

    if (right == -1)
        return -1;
    while (left <= right)
    {
        middle = (left + right)/2;

        if(strcmp(SP[middle].chr,chr)==0)
        {
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
        else if(strcmp(SP[middle].chr,chr)>0) right = middle -1;
        else
            left = middle + 1;
    }
    return left;
}
int exon_cmp(const void *a,const void *b)
{
    struct exon_t *EA,*EB;
    EA = (struct exon_t *)a;
    EB = (struct exon_t *)b;

    if(strcmp(EA->chr,EB->chr)==0)
    {
        if(EA->start==EB->start)
        {
            if (EA->end==EB->end) return 0;
            else if (EA->end<EB->end) return -1;
            else  return 1;
        }
        else if(EA->start < EB->start) return -1;
        else return 1;
    }
    else return strcmp(EA->chr,EB->chr);
}
int find_trans(struct trans_t *SP,int num,char *chr,unsigned int start)
{
    int left = 0, right = num-1, middle;

    if (right == -1)
        return -1;
    while (left <= right)
    {
        middle = (left + right)/2;

        if(strcmp(SP[middle].chr,chr)==0)
        {
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
        else if(strcmp(SP[middle].chr,chr)>0) right = middle -1;
        else
            left = middle + 1;
    }
    return left;
}
int trans_cmp(const void *a,const void *b)
{
    struct trans_t *EA,*EB;
    EA = (struct trans_t *)a;
    EB = (struct trans_t *)b;

    if(strcmp(EA->chr,EB->chr)==0)
    {
        if(EA->start==EB->start)
        {
            if (EA->end==EB->end) return 0;
            else if (EA->end<EB->end) return -1;
            else  return 1;
        }
        else if(EA->start < EB->start) return -1;
        else return 1;
    }
    else return strcmp(EA->chr,EB->chr);
}

struct cigar_t
{
    int l;
    char c;
};
void process_cigar(char *cigarBuf,struct cigar_t *cigar,int *cigar_total,unsigned int *end_pos)
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

                if((cigar[(*cigar_total)].c=='M')||(cigar[(*cigar_total)].c=='X')||(cigar[(*cigar_total)].c=='D')||(cigar[(*cigar_total)].c=='N')) (*end_pos)+= cigar[(*cigar_total)].l;
                (*cigar_total)++;
                j = 0;
            }
        }
    }
}
int read_trans(FILE *gtf,struct trans_t *Trans,struct exon_t *exon)
{
    char f_line[MAX_STRING_LENGTH];
    char name[100];

    int trans_num = 0;
    int exon_num = 0;

    while(fgets(f_line,MAX_STRING_LENGTH,gtf)!=NULL)
    {
        sscanf(f_line,"%s\t%*s\t%*s\t%u\t%u\t%*s\t%*c\t%*c\t%*s %s",exon[exon_num].chr,&(exon[exon_num].start),&(exon[exon_num].end),name);
	    if(trans_num==0)
        {
            strcpy(Trans[trans_num].name,name);
            strcpy(Trans[trans_num].chr,exon[exon_num].chr);
            Trans[trans_num].flag = 0;
            trans_num++;
        }
        else if (strcmp(name,Trans[trans_num-1].name)!=0)
        {
            Trans[trans_num-1].exon_num = exon_num;
            Trans[trans_num-1].exon = (struct exon_t *)malloc(exon_num*sizeof(struct exon_t));

            Trans[trans_num-1].start =  exon[0].start;
            Trans[trans_num-1].end = exon[exon_num-1].end;

            memcpy(Trans[trans_num-1].exon,exon,exon_num*sizeof(struct exon_t));

            strcpy(Trans[trans_num].name,name);
            strcpy(Trans[trans_num].chr,exon[exon_num].chr);
            Trans[trans_num].flag = 0;
            trans_num++;

            memcpy(exon,exon+exon_num,1*sizeof(struct exon_t));
            exon_num=0;
        }

        exon[exon_num].flag = 0;
        exon[exon_num].exact =0;
        exon[exon_num].part = 0;
	    exon_num++;
    }
    Trans[trans_num-1].exon_num = exon_num;
    Trans[trans_num-1].exon = (struct exon_t *)malloc(exon_num*sizeof(struct exon_t));

    Trans[trans_num-1].start =  exon[0].start;
    Trans[trans_num-1].end = exon[exon_num-1].end;

    memcpy(Trans[trans_num-1].exon,exon,exon_num*sizeof(struct exon_t));

    return trans_num;
}
int main(int argc, char *argv[])
{
    FILE *gtf;
    FILE *sam;
    //FILE *write;

    int flag = 0;//0 all 1 best
    if(argv[1][0]=='a') flag = 0;
    else if(argv[1][0]=='b') flag = 1;
    else return 0;

    gtf = fopen(argv[2],"r");
    sam = fopen(argv[3],"r");
    //write = fopen(argv[4],"w");

    char f_line[MAX_STRING_LENGTH];
    char cigar[MAX_STRING_LENGTH];

    unsigned int ref_pos,start_pos,end_pos;
    int read_pos;

    struct cigar_t t_cigar[5000];
    int t_num = 0;
    int i,j,k;

    struct exon_t *exon = (struct exon_t *)malloc(MAX_EXON*sizeof(struct exon_t));
    int exon_num = 0;
    int order = 0;

    struct trans_t *Trans = (struct trans_t *)malloc(MAX_TRANS*sizeof(struct trans_t));
    int trans_num = 0;

    trans_num = read_trans(gtf,Trans,exon);
    qsort(Trans,trans_num,sizeof(struct trans_t),trans_cmp);

    struct exon_t temp_exon[5000];
    int temp_num = 0;

    if(flag==0)
    {
        char chr[400];
        int score = 0,best_score = 0;
        int best_order = 0;
        unsigned int temp_start = 0;

        while(fgets(f_line,MAX_STRING_LENGTH,sam)!=NULL)
        {
            if(f_line[0]=='@') continue;
            sscanf(f_line,"%*s\t%*d\t%s\t%u\t%*d\t%s\t%*s",chr,&ref_pos,cigar);

            process_cigar(cigar,t_cigar,&t_num,&end_pos);
            temp_start=start_pos = ref_pos;

            temp_num = 0;
            for(i = 0;i<t_num;i++)
            {
                switch(t_cigar[i].c)
                {
                    case 'M':
                    case 'X': read_pos+=t_cigar[i].l; ref_pos+=t_cigar[i].l;break;
                    case 'I':
                    case 'S':
                    case 'H': read_pos+=t_cigar[i].l;break;
                    case 'D':
                        if(t_cigar[i].l>20)
                        {
                            temp_exon[temp_num].start = start_pos;
                            temp_exon[temp_num].end = ref_pos;
                            temp_num++;

                            ref_pos+=t_cigar[i].l;
                            start_pos = ref_pos;
                        }
                        else
                            ref_pos+=t_cigar[i].l;
                        break;
                    case 'N':
                        temp_exon[temp_num].start = start_pos;
                        temp_exon[temp_num].end = ref_pos;
                        temp_num++;

                        ref_pos+=t_cigar[i].l;
                        start_pos = ref_pos;
                        break;
                    default:break;
                }
            }
            temp_exon[temp_num].start = start_pos;
            temp_exon[temp_num].end = ref_pos;
            temp_num++;

            best_score = -1;
            best_order = -1;

            order = find_trans(Trans,trans_num,chr,temp_start);
            while(order<trans_num)
            {
                score = 0;

                if(strcmp(chr,Trans[order].chr)!=0) break;
                if(Trans[order].start>temp_exon[temp_num-1].end) break;

                for(j = 0;j<temp_num;j++)
                {
                    for(k =0;k<Trans[order].exon_num;k++)
                    {
                        if(Trans[order].exon[k].start>temp_exon[j].end) break;

                        if((temp_exon[j].start<=Trans[order].exon[k].start)&&(temp_exon[j].end>=Trans[order].exon[k].start))
                            score+= min(temp_exon[j].end-Trans[order].exon[k].start+1,Trans[order].exon[k].end-Trans[order].exon[k].start+1);
                        if((temp_exon[j].end>=Trans[order].exon[k].end)&&(temp_exon[j].start<=Trans[order].exon[k].end))
                            score+= min(Trans[order].exon[k].end-temp_exon[j].start+1,Trans[order].exon[k].end-Trans[order].exon[k].start+1);
                        if((temp_exon[j].start>=Trans[order].exon[k].start)&&(temp_exon[j].end<=Trans[order].exon[k].end))
                            {score+= temp_exon[j].end-temp_exon[j].start+1;break;}
                    }
                }
                if(score>best_score)
                {
                    best_score = score;
                    best_order = order;
                }
                order++;
            }
            if((best_order!=-1)&&(best_score!=0))
            {
                Trans[best_order].flag = 1;
                for(j = 0;j<temp_num;j++)
                {
                    for(k =0;k<Trans[best_order].exon_num;k++)
                    {
                    if(Trans[best_order].exon[k].start>temp_exon[j].end) break;

                    if((temp_exon[j].start<=Trans[best_order].exon[k].start)&&(temp_exon[j].end>=Trans[best_order].exon[k].start)) Trans[best_order].exon[k].part = 1;
                    if((temp_exon[j].end>=Trans[best_order].exon[k].end)&&(temp_exon[j].start<=Trans[best_order].exon[k].end)) Trans[best_order].exon[k].part = 1;
                    if((temp_exon[j].start>=Trans[best_order].exon[k].start)&&(temp_exon[j].end<=Trans[best_order].exon[k].end)) Trans[best_order].exon[k].part = 1;

                    if((temp_exon[j].start<=Trans[best_order].exon[k].start+5)&&(temp_exon[j].start>=Trans[best_order].exon[k].start-5)
                       &&(temp_exon[j].end>=Trans[best_order].exon[k].end-5)&&(temp_exon[j].end<=Trans[best_order].exon[k].end+5)) Trans[best_order].exon[k].exact = 1;
                    }
                }
            }

        }
    }
    else if(flag == 1)
    {
        char temp_name[400];
        char name[400];

        int score = 0,best_score = 0;
        int best_order = 0;
        char best_cigar[MAX_STRING_LENGTH];
        unsigned int best_start = 0;
        char temp_chr[400];
        char chr[400];
        unsigned int temp_start = 0;

        temp_name[0] = '*';
        temp_name[1] = '\0';

        while(fgets(f_line,MAX_STRING_LENGTH,sam)!=NULL)
        {
            if(f_line[0]=='@') continue;
            sscanf(f_line,"%s\t%*d\t%s\t%u\t%*d\t%s\t%*s",name,chr,&temp_start,cigar);

            if(temp_name[0]=='*')
            {
                strcpy(temp_name,name);
                best_score = -1;
                best_order = -1;
                strcpy(temp_chr,chr);
                strcpy(best_cigar,cigar);
                best_start = ref_pos;
            }
            else if(strcmp(temp_name,name)!=0)
            {
                process_cigar(best_cigar,t_cigar,&t_num,&end_pos);
                ref_pos = start_pos = best_start;

                for(i = 0;i<t_num;i++)
                {
                    switch(t_cigar[i].c)
                    {
                        case 'M':
                        case 'X': read_pos+=t_cigar[i].l; ref_pos+=t_cigar[i].l;break;
                        case 'I':
                        case 'S':
                        case 'H': read_pos+=t_cigar[i].l;break;
                        case 'D':
                            if(t_cigar[i].l>20)
                            {
                                temp_exon[temp_num].start = start_pos;
                                temp_exon[temp_num].end = ref_pos;
                                temp_num++;

                                ref_pos+=t_cigar[i].l;
                                start_pos = ref_pos;
                            }
                            else
                                ref_pos+=t_cigar[i].l;
                            break;
                        case 'N':

                            temp_exon[temp_num].start = start_pos;
                            temp_exon[temp_num].end = ref_pos;
                            temp_num++;

                            ref_pos+=t_cigar[i].l;
                            start_pos = ref_pos;
                            break;
                        default:break;
                    }
                }
                temp_exon[temp_num].start = start_pos;
                temp_exon[temp_num].end = ref_pos;
                temp_num++;

                if((best_order!=-1)&&(best_score!=0))
                {
                    Trans[best_order].flag = 1;
                    for(j = 0;j<temp_num;j++)
                    {
                        for(k =0;k<Trans[best_order].exon_num;k++)
                        {
                        if(Trans[best_order].exon[k].start>temp_exon[j].end) break;

                        if((temp_exon[j].start<=Trans[best_order].exon[k].start)&&(temp_exon[j].end>=Trans[best_order].exon[k].start)) Trans[best_order].exon[k].part = 1;
                        if((temp_exon[j].end>=Trans[best_order].exon[k].end)&&(temp_exon[j].start<=Trans[best_order].exon[k].end)) Trans[best_order].exon[k].part = 1;
                        if((temp_exon[j].start>=Trans[best_order].exon[k].start)&&(temp_exon[j].end<=Trans[best_order].exon[k].end)) Trans[best_order].exon[k].part = 1;

                        if((temp_exon[j].start<=Trans[best_order].exon[k].start+5)&&(temp_exon[j].start>=Trans[best_order].exon[k].start-5)
                           &&(temp_exon[j].end>=Trans[best_order].exon[k].end-5)&&(temp_exon[j].end<=Trans[best_order].exon[k].end+5)) Trans[best_order].exon[k].exact = 1;
                        }
                    }
                }
                strcpy(temp_name,name);
                best_score = -1;
                best_order = -1;
                strcpy(temp_chr,chr);
                strcpy(best_cigar,cigar);
                best_start = temp_start;
            }

            process_cigar(cigar,t_cigar,&t_num,&end_pos);
            start_pos = ref_pos = temp_start;

            temp_num = 0;
            for(i = 0;i<t_num;i++)
            {
                switch(t_cigar[i].c)
                {
                    case 'M':
                    case 'X': read_pos+=t_cigar[i].l; ref_pos+=t_cigar[i].l;break;
                    case 'I':
                    case 'S':
                    case 'H': read_pos+=t_cigar[i].l;break;
                    case 'D':
                        if(t_cigar[i].l>20)
                        {
                            temp_exon[temp_num].start = start_pos;
                            temp_exon[temp_num].end = ref_pos;
                            temp_num++;

                            ref_pos+=t_cigar[i].l;
                            start_pos = ref_pos;
                        }
                        else
                            ref_pos+=t_cigar[i].l;
                        break;
                    case 'N':

                        temp_exon[temp_num].start = start_pos;
                        temp_exon[temp_num].end = ref_pos;
                        temp_num++;

                        ref_pos+=t_cigar[i].l;
                        start_pos = ref_pos;
                        break;
                    default:break;
                }
            }
            temp_exon[temp_num].start = start_pos;
            temp_exon[temp_num].end = ref_pos;
            temp_num++;

            order = find_trans(Trans,trans_num,chr,temp_start);
            while(order<trans_num)
            {
                score = 0;

                if(strcmp(chr,Trans[order].chr)!=0) break;
                if(Trans[order].start>temp_exon[temp_num-1].end) break;

                for(j = 0;j<temp_num;j++)
                {
                    for(k =0;k<Trans[order].exon_num;k++)
                    {
                        if(Trans[order].exon[k].start>temp_exon[j].end) break;

                        if((temp_exon[j].start<=Trans[order].exon[k].start)&&(temp_exon[j].end>=Trans[order].exon[k].start))
                            score+= min(temp_exon[j].end-Trans[order].exon[k].start+1,Trans[order].exon[k].end-Trans[order].exon[k].start+1);
                        if((temp_exon[j].end>=Trans[order].exon[k].end)&&(temp_exon[j].start<=Trans[order].exon[k].end))
                            score+= min(Trans[order].exon[k].end-temp_exon[j].start+1,Trans[order].exon[k].end-Trans[order].exon[k].start+1);
                        if((temp_exon[j].start>=Trans[order].exon[k].start)&&(temp_exon[j].end<=Trans[order].exon[k].end))
                            {score+= temp_exon[j].end-temp_exon[j].start+1;break;}
                    }
                }
                if(score>best_score)
                {
                    best_score = score;
                    strcpy(temp_chr,chr);
                    strcpy(best_cigar,cigar);
                    best_start = temp_start;
                    best_order = order;
                }
                order++;
            }
        }
    }

    int find_trans = 0;
    exon_num = 0;

    for(i = 0;i<trans_num;i++)
    {
        if(Trans[i].flag==1) find_trans++;

        memcpy(exon+exon_num,Trans[i].exon,Trans[i].exon_num*sizeof(struct exon_t));
        exon_num+=Trans[i].exon_num;
    }
    qsort(exon,exon_num,sizeof(struct exon_t),exon_cmp);

    int start = 0;
    exon[start].flag = 1;
    for(i = 1;i<exon_num;i++)
    {
        if((exon[i].start==exon[start].start)&&(exon[i].end==exon[start].end)&&(strcmp(exon[i].chr,exon[start].chr)==0))
        {
            exon[i].flag = 0;
            exon[start].part|=exon[i].part;
            exon[start].exact|=exon[i].exact;
        }
        else{start = i;exon[start].flag = 1;}
    }

    int exon15 = 0;
    int exon20 = 0;
    int exon30 = 0;
    int exon40 = 0;
    int exon50 = 0;
    int exon60 = 0;
    int exon70 = 0;
    int exon80 = 0;
    int exon90 = 0;
    int exon100 = 0;
    int length = 0;
    int total = 0;
    int find = 0;
    for(i = 0;i<exon_num;i++)
    {
        if(exon[i].flag==0) continue;
        total++;

        length = exon[i].end-exon[i].start+1;
        if(exon[i].part==0)
        {
            //fprintf(write,"%s\t%u\t%u\t%d\n",exon[i].chr,exon[i].start,exon[i].end,length);
            continue;
        }
	find ++;
        if(length<16) exon15++;
        else if(length<21) exon20++;
        else if(length<31) exon30++;
        else if(length<41) exon40++;
        else if(length<51) exon50++;
        else if(length<61) exon60++;
        else if(length<71) exon70++;
        else if(length<81) exon80++;
        else if(length<91) exon90++;
        else exon100++;
    }
    printf("total trans:   %d\n",trans_num);
    printf("find trans:     %d\n",find_trans);
    printf("total exons	%d\n",total);
    printf("find exons	%d\n",find);
    printf("find exons length <=15	%d\n",exon15);
    printf("find exons length <=20&>15	%d\n",exon20);
    printf("find exons length <=30&>20	%d\n",exon30);
    printf("find exons length <=40&>30	%d\n",exon40);
    printf("find exons length <=50&>40	%d\n",exon50);
    printf("find exons length <=60&>50	%d\n",exon60);
    printf("find exons length <=70&>60	%d\n",exon70);
    printf("find exons length <=80&>70	%d\n",exon80);
    printf("find exons length <=90&>80	%d\n",exon90);
    printf("find exons length >90	%d\n",exon100);


    free(exon);

    fclose(gtf);
    fclose(sam);
    //fclose(write);

    return 0;
}
