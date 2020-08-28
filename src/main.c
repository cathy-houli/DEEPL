#include "main.h"

int usage()
{
    fprintf(stdout, "\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "Program: use RNA-seq to find junctions, use align file to find more junctions\n");
    fprintf(stdout, "Version: %s\n", PACKAGE_VERSION);
    fprintf(stdout, "Contact: Hou.L <ye.wenzhu@gmail.com>\n\n");

    fprintf(stdout, "Usage:   DEEP <command>\n\n");

    fprintf(stdout, "command :\n");
    fprintf(stdout, "   index:         build reference hash index\n");
    fprintf(stdout, "   micro:         run micro mode\n");
    fprintf(stdout, "   complete:      run complete mode\n");
    fprintf(stdout, "   version:       show current version\n\n");

    fprintf(stdout, "index:  DEEP index input_path output_path\n\n");

    fprintf(stdout, "   input_path:    reference files path, each need file's name must be *.fa\n");
    fprintf(stdout, "   output_path:   need enough space\n\n");


    fprintf(stdout, "complete option :\n");
    fprintf(stdout, "<basic> :\n");
    fprintf(stdout, "   -B <STRING>:   reference BWT index path\n");
    fprintf(stdout, "   -H <STRING>:   hash files path, each need file's name must be *.hash with *.ann files\n");
    fprintf(stdout, "   -f         :   fa formate\n");
    fprintf(stdout, "   -q         :   fq formate\n");
    fprintf(stdout, "   -1 <STRING>:   input files, need *.fa,*.fq files,split with ','\n");
    fprintf(stdout, "   -2 <STRING>:   input pair-end files, need *.fa,*.fq files,split with ','\n");
    fprintf(stdout, "   -O <STRING>:   output path,need enough space\n\n");

    fprintf(stdout, "<alignment parameter> :\n");
    fprintf(stdout, "   -t <INT>:      output read score threshold(90)\n");
    fprintf(stdout, "   -r <INT>:      search area(500000)\n");
    fprintf(stdout, "   -a :           print all read sites while coverage score and alignment score over set threld will be used(-b)\n");
    fprintf(stdout, "   -b :           print best read sites while coverage score and alignment score over set threld will be used(-b)\n");
    fprintf(stdout, "   -p <INT>:      thread(1)\n");
    fprintf(stdout, "   -m <INT>:      match score(1)\n");
    fprintf(stdout, "   -s <INT>:      miss score(-1)\n");
    fprintf(stdout, "   -g <INT>:      gap score(-3)\n");
    fprintf(stdout, "   -e <INT>:      min exon length(12)\n");
    fprintf(stdout, "   -d             deep mode(0)\n\n");

    return 1;
}
char base2char[5] = {'A','C','G','T'};
unsigned char nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
void opt_init(struct m_opt **opt)
{
    (*opt)->idx = (bwaidx_t*)calloc(1, sizeof(bwaidx_t));

    (*opt)->input_file_1 = (struct file_list *)calloc(1,sizeof(struct file_list));
	(*opt)->input_file_1->file = (struct file_name *)calloc(MAX_FILE_NUM,sizeof(struct file_name));

	(*opt)->file_flag = 0;

	(*opt)->chr = (struct chr_list *)calloc(1,sizeof(struct chr_list));
	(*opt)->chr->list = (struct chr_t *)calloc(MAX_FILE_NUM,sizeof(struct chr_t));

	(*opt)->SP_MAX = 100000;
	(*opt)->SP_block_MAX = 100000;
	(*opt)->SP_num = 0;
	(*opt)->SP = (struct splice_list *)calloc((*opt)->SP_MAX,sizeof(struct splice_list));

	(*opt)->change_length = 5;

	(*opt)->input_mode = FQ_FILE;
	(*opt)->area = 1000000;

	(*opt)->match = 2;
	(*opt)->miss = 4;
	(*opt)->insert = 2;
	(*opt)->gap = 2;
	(*opt)->splice = 20;

	(*opt)->score_t= 0.9;

	(*opt)->thread_num = 1;
	(*opt)->deep_mode = 0;
	(*opt)->pass = 1;

	(*opt)->cut_tail = 0;

	(*opt)->thread_block = READ_BUF_LENGTH;
	(*opt)->result_block = (*opt)->thread_block*5;

	(*opt)->un_num = 0;
}
void opt_free(struct m_opt **opt)
{
    free((*opt)->idx);

    free((*opt)->input_file_1->file);
    free((*opt)->input_file_1);

    unsigned int i = 0;
    for(i = 0;i<(*opt)->chr->total;i++)
        free((*opt)->chr->list[i].seq);
    free((*opt)->chr->list);
	free((*opt)->chr);

	for (i=0; i<pow(4,KMER_FRONT); i++)
    {
        if((*opt)->c_hash[i]!=NULL) free((*opt)->c_hash[i]);
    }
    free((*opt)->c_num);

	free((*opt)->SP);
}
int check_inputfile(struct file_list *input_list,char *input_file_name)
{
    char *p;
    int i = 0;
    p=strtok(input_file_name,",");
    while(p)
    {
        if((access(p,0))== -1)
        {
            fprintf(stdout,"error!!! cannot find file %s...\n",p);
            return 1;
        }

        if(i<MAX_FILE_NUM)
        {
            strncpy(input_list->file[i].name,p,MAX_NAME_LENGTH);
            i++;
        }
        else
        {
            fprintf(stdout,"error!!! too many input files...\n");
            return 1;
        }
        p = strtok(NULL,",");
    }
    input_list->total = i;
    return 0;
}
int check_hash_index(struct m_opt *opt)
{
    char name[MAX_STRING_LENGTH];
    sprintf(name,"%s/chr.hash",opt->Hash_path);
    if((access(name,0))== -1)
    {
        fprintf(stdout,"error!!! cannot find %s...\n",name);
        return 1;
    }
    return 0;
}
int main(int argc, char *argv[])
{
    if(argc<3) return usage();

    time_t present_time;
    struct tm *present_tm;

    time(&present_time);
    present_tm = localtime(&present_time);
    fprintf(stdout,"[DEEP-LONG]begin micro program...%s\n",asctime(present_tm));

    if(strcmp(argv[1],"index")==0) return index_main(argc-1, argv+1);
    else if(strcmp(argv[1],"aln")!=0) return usage();

    opt = (struct m_opt *)calloc(1,sizeof(struct m_opt));
    opt_init(&opt);

    char input_file_name[MAX_STRING_LENGTH];
    input_file_name[0]='\0';

    int i = 0;
    for(i = 2;i<argc;i++)
    {
        if(argv[i][0]=='-')
        {
        switch (argv[i][1]) {
            case 'B': strncpy(opt->BWTpath,argv[++i],MAX_NAME_LENGTH);break;
            case 'H': strncpy(opt->Hash_path,argv[++i],MAX_NAME_LENGTH);break;
            case 'O': strncpy(opt->Output_path,argv[++i],MAX_NAME_LENGTH);break;
            case '1': strncpy(input_file_name,argv[++i],MAX_NAME_LENGTH);break;
            case 'r': opt->area = atoi(argv[++i]); break;
            case 's': opt->miss = atoi(argv[++i]); break;
            case 'm': opt->match = atoi(argv[++i]); break;
            case 'g': opt->gap = atoi(argv[++i]); break;
            case 'p': opt->thread_num = atoi(argv[++i]); break;
            case 'f': opt->input_mode = FA_FILE; break;
            case 'q': opt->input_mode = FQ_FILE; break;
            case 'd': opt->deep_mode = 1; break;
            case 't': opt->cut_tail = atoi(argv[++i]); break;
            default: return usage();
            }
        }
		else return usage();
    }
    opt->pass = 1;
    opt->thread_block = min(1500,READ_BUF_LENGTH*opt->thread_num);

    if(input_file_name[0]=='\0') return usage();
    else {if(check_inputfile(opt->input_file_1,input_file_name)) return 1;}

    //if(check_outpath(opt)) return 1;
    if(check_hash_index(opt)) return 1;

    time(&present_time);
    present_tm = localtime(&present_time);
    fprintf(stdout,"[DEEP-LONG]load index...%s\n",asctime(present_tm));
    if(load_index(opt)) return 1;
    if(load_hash(opt)) return 1;

    opt->l_pac = opt->idx->bns->l_pac;

    time(&present_time);
    present_tm = localtime(&present_time);
    fprintf(stdout,"[DEEP-LONG]begin align...%s\n",asctime(present_tm));

    if(seed_align(opt)) return 1;

    time(&present_time);
    present_tm = localtime(&present_time);
    fprintf(stdout,"[DEEP-LONG]recheck splice...%s\n",asctime(present_tm));
    if(re_check(opt)) return 1;

    bwa_idx_destroy(opt->idx);

    opt_free(&opt);
    free(opt);

    time(&present_time);
    present_tm = localtime(&present_time);
    fprintf(stdout,"[DEEP-LONG]exit program...%s\n",asctime(present_tm));
    return 0;
}

