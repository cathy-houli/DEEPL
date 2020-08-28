#include "main.h"
KSEQ_INIT(gzFile, gzread)

#define MAX_CHR_HASH_SAME 2*1024*1024*256//500
#define BUCKET_LENGTH (100<<20)
#define MAX_INDEX_LENGTH 3LL*1300*1024*1024
#define CHR_MARK ((1LL<<((KMER_FRONT)<<1))-1)
struct heap_array
{
    unsigned int *heap;
    unsigned int *heap_num;
    unsigned int *order;
    unsigned int heap_length;
};
struct heap_array *heap;
#define BUCKET_NUM 100
struct sort_chr
{
    uint16_t mark;
    unsigned int site;
};
int chr_hash_cmp(const void *a,const void *b)
{
    struct sort_chr *EA,*EB;
    EA = (struct sort_chr *)a;
    EB = (struct sort_chr *)b;

    if(EA->mark==EB->mark)
    {
        if(EA->site==EB->site) return 0;
        if(EA->site<EB->site) return -1;
        if(EA->site>EB->site) return 1;
    }
    else if(EA->mark<EB->mark) return -1;
    else if(EA->mark>EB->mark) return 1;
    return 0;
}
struct chr_sort_hash
{
    struct sort_chr *buf;
    unsigned int Enum;
    unsigned int Bnum;
};
struct chr_sort_hash *sort_array;
struct sort_chr *chr_buf;
void SwapHeap(unsigned int a,unsigned int b,unsigned int *heap,unsigned int *heap_num,unsigned int *order)
{
    unsigned int temp = 0;
    temp = heap[a];
    heap[a] = heap[b];
    heap[b] = temp;

    temp = heap_num[a];
    heap_num[a] = heap_num[b];
    heap_num[b] = temp;

    temp = order[a];
    order[a] = order[b];
    order[b] = temp;
}
void HeapAdjust_Chr(unsigned int *heap,unsigned int *heap_num,unsigned int *order,unsigned int heap_length,unsigned int i)  //调整堆
{
    int LC=2*i;       //i的左孩子节点序号
    int RC=2*i+1;     //i的右孩子节点序号
    int min=i;            //临时变量
    if(i<=heap_length/2)          //如果i是叶节点就不用进行调整
    {
        if((LC<=heap_length)
        &&((chr_buf[heap[LC]].mark<chr_buf[heap[min]].mark)||(chr_buf[heap[LC]].mark==chr_buf[heap[min]].mark&&chr_buf[heap[LC]].site<chr_buf[heap[min]].site)))
            min=LC;
        if((RC<=heap_length)
        &&((chr_buf[heap[RC]].mark<chr_buf[heap[min]].mark)||(chr_buf[heap[RC]].mark==chr_buf[heap[min]].mark&&chr_buf[heap[RC]].site<chr_buf[heap[min]].site)))
            min=RC;
        if(min!=i)
        {
            SwapHeap(i,min,heap,heap_num,order);
            HeapAdjust_Chr(heap,heap_num,order,heap_length,min); //避免调整之后以max为父节点的子树不是堆
        }
    }
}
void ReHeap_Chr(unsigned int *heap,unsigned int *heap_num,unsigned int *order,unsigned int *heap_length,FILE *InPut_file[])
{
    heap_num[1]--;
    if(heap_num[1]>0)
    {
        if(fread(&chr_buf[order[1]],sizeof(struct sort_chr),1,InPut_file[order[1]])!=1)
        {
            heap_num[1] = 0;
            SwapHeap(1,*heap_length,heap,heap_num,order);
            (*heap_length)--;
        }
    }
    else{
        SwapHeap(1,*heap_length,heap,heap_num,order);
        (*heap_length)--;
    }
    HeapAdjust_Chr(heap,heap_num,order,*heap_length,1);
}
int index_main(int argc, char *argv[])
{
    if(argc<3)return usage();

    time_t present_time;
    struct tm *present_tm;

    time(&present_time);
    present_tm = localtime(&present_time);
    fprintf(stdout,"build hash index...%s\n",asctime(present_tm));

    char system_order[MAX_STRING_LENGTH];

    if((access(argv[1],0))== -1)
    {
        fprintf(stdout,"cannot access %s, please check your input.\n",argv[1]);
        return 1;
    }
    if((access(argv[2],0))== -1)
    {
        fprintf(stdout,"cannot access %s, create output files folder %s...\n",argv[2],argv[2]);
        sprintf(system_order,"mkdir %s",argv[2]);
        system(system_order);

        if((access(argv[2],0))== -1)
        {
            fprintf(stdout,"error!!! cannot create output files folder...\n");
            return 1;
        }
    }

    strcpy(system_order,"ls ");
    strcat(system_order,argv[1]);
    strcat(system_order,"/*.fa >");
    strcat(system_order,argv[2]);
    strcat(system_order,"/ref_name.txt");
    system(system_order);

    FILE *Chr_info;
    FILE *Name_info;
    FILE *Output_hash;

    gzFile ref;
    kseq_t *q;

    char chr_file[MAX_NAME_LENGTH];
    char name_file[MAX_NAME_LENGTH];
    char in_file[MAX_NAME_LENGTH];
    char out_file[MAX_NAME_LENGTH];

    strcpy(chr_file,argv[2]);
    strcat(chr_file,"/chr_list.txt");

    strcpy(name_file,argv[2]);
    strcat(name_file,"/ref_name.txt");

    Chr_info = fopen(chr_file,"w");
    if (Chr_info == NULL)
        return 1;

    Name_info = fopen(name_file,"r");
    if (Name_info == NULL)
        return 1;

    strcpy(system_order,"sort -k1,1 ");
    strcat(system_order,name_file);
    strcat(system_order," -o ");
    strcat(system_order,name_file);
    system(system_order);

    char f_line[MAX_STRING_LENGTH];
    uint32_t base = 0;
    uint32_t temp_base = -1;
    unsigned int i = 0;
    unsigned int rel_site = 0;
    unsigned int read_site = 0;
    int j = 0,move = 0;
    int N = 0;

    unsigned int *buffer=(unsigned int *)calloc(MAX_CHR_HASH_SAME,sizeof(unsigned int));
    int buffer_num = 0;
    int buffer_flag = 0;

    int *num;
    num= (int *)calloc(pow(4,KMER_FRONT),sizeof(int));

    chr_buf = (struct sort_chr *)calloc(MAX_INDEX_LENGTH/BUCKET_LENGTH+1,sizeof(struct sort_chr));
    sort_array = (struct chr_sort_hash *)calloc(1,sizeof(struct chr_sort_hash));
    sort_array->buf = (struct sort_chr *)calloc(BUCKET_LENGTH,sizeof(struct sort_chr));
    sort_array->Enum = 0;
    sort_array->Bnum = 0;

    heap = (struct heap_array *)calloc(1,sizeof(struct heap_array));
    heap->heap = (unsigned int *)calloc(MAX_INDEX_LENGTH/BUCKET_LENGTH+1,sizeof(unsigned int));
    heap->heap_num = (unsigned int *)calloc(MAX_INDEX_LENGTH/BUCKET_LENGTH+1,sizeof(unsigned int));
    heap->order = (unsigned int *)calloc(MAX_INDEX_LENGTH/BUCKET_LENGTH+1,sizeof(unsigned int));
    heap->heap_length = 0;

    FILE *temp_file;
    char temp_name[MAX_NAME_LENGTH];

    int heap_length = 0;
    int l = 0;

    while (fgets(f_line,MAX_STRING_LENGTH,Name_info)!=NULL)
    {
        sscanf(f_line,"%s",in_file);
        ref = gzopen(in_file,"r");
        q = kseq_init(ref);

        while ((l = kseq_read(q)) > 0)
        {
        time(&present_time);
        present_tm = localtime(&present_time);

        fprintf(stdout,"process %s...%s",q->name.s,asctime(present_tm));
        fprintf(Chr_info,"%s\t%u\t%u\n",q->name.s,(unsigned int)q->seq.l,rel_site);

        move = 0;
        base = 0;

        for (i = 0; i < q->seq.l-KMER_FRONT; i++)
        {
            N = 0;
            if(sort_array->Enum>0)
            {
                move =i+rel_site-sort_array->buf[sort_array->Enum-1].site;
                base = sort_array->buf[sort_array->Enum-1].mark;
            }
            if((move>0)&&(move<(KMER_FRONT))) j = KMER_FRONT-move;
            else j = 0;
            for (; j <KMER_FRONT; j++)
            {
                if(nst_nt4_table[(int)q->seq.s[i+j]] >3)
                {
                    N = 1;
                    break;
                }
                base = base << 2 | nst_nt4_table[(int)q->seq.s[i+j]];
            }
            if(N) continue;
            base = base & CHR_MARK;

            if(i%2==0)
            {
                sort_array->buf[sort_array->Enum].mark = base;
                sort_array->buf[sort_array->Enum].site = i+rel_site;
                sort_array->Enum++;

                if(sort_array->Enum>=BUCKET_LENGTH)
                {
                qsort(sort_array->buf,BUCKET_LENGTH,sizeof(struct sort_chr),chr_hash_cmp);

                sprintf(temp_name,"%s/%d.temp",argv[2],sort_array->Bnum);
                temp_file = fopen(temp_name,"wb");
                if(temp_file == NULL)
                    fprintf(stdout,"cannot open file %s...\n",temp_name);

                if(fwrite(sort_array->buf,sizeof(struct sort_chr),BUCKET_LENGTH,temp_file)!=BUCKET_LENGTH)
                    fprintf(stdout,"cannot write file...\n");
                fclose(temp_file);
                sort_array->Enum = 0;
                sort_array->Bnum++;
                }
            }
        }
        rel_site += q->seq.l;
        }

        gzclose(ref);
    }
    if(sort_array->Enum!=0)
    {
        qsort(sort_array->buf,sort_array->Enum,sizeof(struct sort_chr),chr_hash_cmp);

        sprintf(temp_name,"%s/%d.temp",argv[2],sort_array->Bnum);
        temp_file = fopen(temp_name,"wb");
        if(temp_file == NULL)
            fprintf(stdout,"cannot open file %s...\n",temp_name);

        if(fwrite(sort_array->buf,sizeof(struct sort_chr),sort_array->Enum,temp_file)!=sort_array->Enum)
            fprintf(stdout,"cannot write file...\n");
        fclose(temp_file);
    }
    free(sort_array->buf);

    sprintf(out_file,"%s/chr.hash",argv[2]);
    Output_hash = fopen(out_file,"w");

    FILE *InPut_file[MAX_INDEX_LENGTH/BUCKET_LENGTH+1];

    //float score = 0;
    for(i = 0;i<sort_array->Bnum;i++)
    {
        sprintf(in_file,"%s/%d.temp",argv[2],i);
        InPut_file[i] = fopen(in_file,"rb");
        if (InPut_file[i] == NULL)
        {
            fprintf(stdout,"cannot open file %s...\n",in_file);
            continue;
        }
        if(fread(&chr_buf[i],sizeof(struct sort_chr),1,InPut_file[i])!=1)
        {
            fprintf(stdout,"cannot read seq...\n");
            continue;
        }

        heap->heap[i+1] = i;
        heap->heap_num[i+1] = BUCKET_LENGTH;
        heap->order[i+1] = i;
    }
    heap_length =heap->heap_length = sort_array->Bnum;
    if(sort_array->Enum!=0)
    {
        sprintf(in_file,"%s/%d.temp",argv[2],sort_array->Bnum);
        InPut_file[sort_array->Bnum] = fopen(in_file,"rb");
        if (InPut_file[sort_array->Bnum] == NULL)
        {
            fprintf(stdout,"cannot open file %s...\n",in_file);
        }
        if(fread(&chr_buf[sort_array->Bnum],sizeof(struct sort_chr),1,InPut_file[sort_array->Bnum])!=1)
        {
            fprintf(stdout,"cannot read seq...\n");
        }

        heap->heap[sort_array->Bnum+1] = sort_array->Bnum;
        heap->heap_num[sort_array->Bnum+1] = sort_array->Enum;
        heap->order[sort_array->Bnum+1] = sort_array->Bnum;
        heap_length =heap->heap_length = sort_array->Bnum+1;
    }
    for(i=heap->heap_length/2;i>=1;i--)
        HeapAdjust_Chr(heap->heap,heap->heap_num,heap->order,heap->heap_length,i);

    int r = heap->heap[1];

    memset(num,0,sizeof(int)*pow(4,KMER_FRONT));

    ReHeap_Chr(heap->heap,heap->heap_num,heap->order,&heap->heap_length,InPut_file);

    int x = 0;
    temp_base = -1;
    while (heap->heap_num[1]>0)
    {
        r = heap->heap[1];
        base = chr_buf[r].mark;
        read_site = chr_buf[r].site;

        if (base!=temp_base)
        {
            if (temp_base == -1)
            {
                    temp_base = base;

                    buffer_num = 0;
                    buffer_flag = 0;

                    buffer[buffer_num] = read_site;
                    buffer_num++;
            }
            else
            {
                    if (buffer_flag == 0)
                    {
                        num[temp_base] = buffer_num;
                        fprintf(Output_hash,">%u %d\n",temp_base,buffer_num);
                        for(x= 0;x<buffer_num;x++)
                            fprintf(Output_hash,"%u\n",buffer[x]);
                    }
                    buffer_num = 0;
                    buffer_flag = 0;

                    temp_base = base;

                    buffer[buffer_num] = read_site;
                    buffer_num++;
            }
        }
        else
        {
                if (buffer_flag == 1)
                {
                    ReHeap_Chr(heap->heap,heap->heap_num,heap->order,&heap->heap_length,InPut_file);
                    continue;
                }
                if ((buffer_num+1)>MAX_CHR_HASH_SAME)
                {
                    buffer_flag = 1;
                    buffer_num = 0;
                }
                else
                {
                    buffer[buffer_num] = read_site;
                    buffer_num++;
                }
        }
        ReHeap_Chr(heap->heap,heap->heap_num,heap->order,&heap->heap_length,InPut_file);
    }
    if (buffer_flag == 0)
    {
        num[temp_base] = buffer_num;
        fprintf(Output_hash,">%u %d\n",temp_base,buffer_num);
        for(x= 0;x<buffer_num;x++)
            fprintf(Output_hash,"%u\n",buffer[x]);
    }
    fclose(Output_hash);

    for(i = 0;i<heap_length;i++)
    {
        fclose(InPut_file[i]);
        sprintf(in_file,"%s/%d.temp",argv[2],i);
        remove(in_file);
    }

    fclose(Chr_info);
    fclose(Name_info);

    remove(name_file);
    kseq_destroy(q);
    free(num);

    free(heap->heap);
    free(heap->heap_num);
    free(heap->order);
    free(heap);
    free(chr_buf);

    free(buffer);
    free(sort_array);

    time(&present_time);
    present_tm = localtime(&present_time);
    fprintf(stdout,"complete build hash index...%s\n",asctime(present_tm));
    return 0;
}
uint16_t char2base_16(char *readseq,int length)
{
    uint16_t baseseq = 0;
    int i, base = 0;

    for (i = 0; i < length; i++)
    {
        base = nst_nt4_table[(int)readseq[i]];
        if(base == -1)
            return -1;
        baseseq = baseseq << 2 | base;
    }
    return baseseq;
}
int hash_c_find(unsigned int **hash_front,int *hash_num,uint16_t base,unsigned int site)
{
    int left = 0, right = hash_num[base]-1, middle;

    if (right == -1)
        return 0;
    while (left <= right)
    {
        middle = (left + right)/2;

        if (hash_front[base][middle] == site)
        {
            while ((middle-1>=0)&&(hash_front[base][middle-1] == site))
                middle--;
            return middle;
        }
        else if (hash_front[base][middle] > site)
            right = middle -1;
        else
            left = middle + 1;
    }
    return left;
}
int insert_c_hash(unsigned int **hash_front,int *hash_num,uint16_t base,FILE *hash_file)
{
    unsigned int *temp;
    temp = (unsigned int *)calloc(hash_num[base], sizeof(unsigned int));//申请空间
    if(temp == NULL)
    {
        fprintf(stdout,"There have no enough space...\n");
        return 1;
    }

    int i;
    char f_line[MAX_STRING_LENGTH];
    for (i = 0;i<hash_num[base];i++)
    {
        if(fgets(f_line,MAX_STRING_LENGTH,hash_file)!=NULL)
            sscanf(f_line,"%u",temp+i);
        else return 1;
    }

    hash_front[(int)base] = temp;
    return 0;
}
int build_c_hash(struct m_opt *opt,FILE *hash_file)
{
    uint32_t i = 0;
    int num = 0;
    char f_line[MAX_STRING_LENGTH];
    while (fgets(f_line,MAX_STRING_LENGTH,hash_file)!=NULL)
    {
        if(f_line[0]=='>')
        {
            sscanf(f_line,">%u %d",&i,&num);
            opt->c_num[i]=num;
        }
        if(insert_c_hash(opt->c_hash,opt->c_num,i,hash_file)) return 1;
    }
    return 0;
}
int load_hash(struct m_opt *opt)
{
    FILE *Hash_file;
    char name[MAX_STRING_LENGTH];
    sprintf(name,"%s/chr.hash",opt->Hash_path);
    Hash_file = fopen(name,"r");

    opt->c_hash = (unsigned int **)calloc(pow(4,KMER_FRONT),sizeof(unsigned int *));
	opt->c_num= (int *)calloc(pow(4,KMER_FRONT),sizeof(int));

    if(build_c_hash(opt,Hash_file)) return 1;
    fclose(Hash_file);
    return 0;
}

