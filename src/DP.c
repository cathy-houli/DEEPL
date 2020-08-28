#include "main.h"

int Tail_DP(char *chr,unsigned int start,unsigned int end,int mode,char* string,int length,struct cigar_t *cigar,int *cigar_num,struct Splice_DP_t *DP,unsigned int *start_pos,int S)//0 to front 1 to back
{
    if(length==0)
    {
        (*cigar_num) = 0;return 0;
    }
    //if(length>100) {(*cigar_num) = 0;return 0;}
    int GAP = 4;
    int MaxPenalty = -65536;

    char *text = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
    char *ref = (char *)calloc(MAX_READ_LENGTH, sizeof(char));

    struct DPC_t *cigarD = DP->cigarD;
    int DN = 0;

    int i, j,x;

    int  n = max(length*(1.1),length+10);
    if (mode==0) n = min((start-end),n);
    else n = min((end-start),n);

    int  m = length;

    if((n>=MAX_DP_AREA)||(m>=MAX_DP_AREA))
    {
        if(S)
        {
            (*cigar_num) = 1;
            cigar[0].c = 'S';
            cigar[0].l = length;
        }
        else
            (*cigar_num) = 0;
        free(ref);
        free(text);
        return -10000;
    }

    if(mode==0)
    {
        for(i = 0;i<n;i++)
            ref[i] = chr[start-i];
        ref[n] = '\n';

        for(i = 0;i<m;i++)
            text[i] = string[m-i-1];
        text[m] = '\0';
    }
    else
    {
        for(i = 0;i<n;i++)
            ref[i] = chr[start+i];
        ref[n] = '\n';

        for(i = 0;i<m;i++)
            text[i] = string[i];
        text[m] = '\0';
    }


	int** r =  DP->r;
	int** t =  DP->t;
	int** d =  DP->d;

	// initialization
	r[0][0] = t[0][0] = d[0][0] = 0;
	for (i = 1; i <= m; i++)
	{
		r[i][0] = MaxPenalty;
		d[i][0] = t[i][0] = -(GAP + (i-1)*opt->gap);
	}

	for (j = 1; j <= n; j++)
	{
		t[0][j] = MaxPenalty;
		d[0][j] = r[0][j] =  -min(opt->splice,(GAP + (j-1)*opt->gap));
	}

	for (i = 1; i <= m; i++)
	{
		for (j = 1; j <= n; j++)
		{
			r[i][j] = max(r[i][j - 1] - opt->gap, d[i][j - 1] - GAP);
			t[i][j] = max(t[i - 1][j] - opt->gap, d[i - 1][j] - GAP);
			//if (s1[i - 1] == 'N' || s2[j - 1] == 'N') s[i][j] = max(s[i - 1][j - 1] + 1, r[i][j], t[i][j]);
			//else s[i][j] = max(s[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 1 : -1), r[i][j], t[i][j]);
			d[i][j] = max(d[i - 1][j - 1] + (nst_nt4_table[(int)text[i - 1]] == nst_nt4_table[(int)ref[j - 1]] ? opt->match : -opt->miss), max(r[i][j], t[i][j]));
		}
	}

	int y = 0;
    for(j = 1;j<=n;j++)
    {
        if(d[m][j]>d[m][y]) y = j;
    }
	i = m;j = y;
	DN = 0;
	int max_x = 0;

	while (i > 0 || j > 0)
    {
		if (d[i][j] == r[i][j])
        {
            cigarD[DN].cigar = 'D';
            cigarD[DN].score = d[i][j];
            cigarD[DN].m = i;
            cigarD[DN].n = (mode==0)?(start-j+1):(j+start-1);
            if(cigarD[DN].score>cigarD[max_x].score)
                max_x = DN;
            DN++;
			j--;
		}
		else if (d[i][j] == t[i][j])
        {
			cigarD[DN].cigar = 'I';
            cigarD[DN].score = d[i][j];
            cigarD[DN].m = i;
            cigarD[DN].n = (mode==0)?(start-j+1):(j+start-1);
            if(cigarD[DN].score>cigarD[max_x].score)
                max_x = DN;
            DN++;
			i--;
		}
		else
        {
            if(d[i][j]>d[i-1][j-1]) cigarD[DN].cigar = 'M';
            else cigarD[DN].cigar = 'X';
            cigarD[DN].score = d[i][j];
            cigarD[DN].m = i;
            cigarD[DN].n = (mode==0)?(start-j+1):(j+start-1);
            if(cigarD[DN].score>cigarD[max_x].score)
                max_x = DN;
            DN++;
			i--, j--;
		}
	}

	if(!S) max_x = 0;

    x = 0;
    cigar[x].c = cigarD[DN-1].cigar;
    cigar[x].l = 0;
    for(i=DN-1;i>=max_x;i--)
    {
        if(cigar[x].c==cigarD[i].cigar)
            cigar[x].l++;
        else
        {
            x++;
            cigar[x].c = cigarD[i].cigar;
            cigar[x].l = 1;
        }
    }
    if(cigar[x].l != 0)x++;

    if(cigarD[max_x].m<length)
    {
        cigar[x].c = 'S';
        cigar[x].l = length-cigarD[max_x].m;
        x++;
    }
    (*cigar_num) = x;

    *start_pos = cigarD[max_x].n;

    char c;
    int l;
    if(mode==0)
    {
        for(i = 0;i<x/2;i++)
        {
            c = cigar[i].c;
            l = cigar[i].l;

            cigar[i].c = cigar[x-i-1].c;
            cigar[i].l = cigar[x-i-1].l;

            cigar[x-i-1].c = c;
            cigar[x-i-1].l = l;
        }
    }
    free(text);
    free(ref);

	return cigarD[max_x].score;
}
int NW(char *ref,int n,char* text,int m,struct cigar_t *cigar,int *cigar_num,struct Splice_DP_t *DP)
{
    int GAP = 4;
    int MaxPenalty = -65536;

    int i, j,x;

	int** r =  DP->r;
	int** t =  DP->t;
	int** d =  DP->d;

	if((n>=MAX_DP_AREA)||(m>=MAX_DP_AREA))
    {
        (*cigar_num) = 0;
        return -10000;
    }

	// initialization
	r[0][0] = t[0][0] = d[0][0] = 0;
	for (i = 1; i <= m; i++)
	{
		r[i][0] = MaxPenalty;
		d[i][0] = t[i][0] = -min(65536,(GAP + (i-1)*opt->gap));
	}

	for (j = 1; j <= n; j++)
	{
		t[0][j] = MaxPenalty;
		d[0][j] = r[0][j] = -min(65536,(GAP + (j-1)*opt->gap));
	}

	for (i = 1; i <= m; i++)
	{
		for (j = 1; j <= n; j++)
		{
			r[i][j] = max(r[i][j - 1] - opt->gap, d[i][j - 1] - GAP);
			t[i][j] = max(t[i - 1][j] - opt->gap, d[i - 1][j] - GAP);
			//if (s1[i - 1] == 'N' || s2[j - 1] == 'N') s[i][j] = max(s[i - 1][j - 1] + 1, r[i][j], t[i][j]);
			//else s[i][j] = max(s[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 1 : -1), r[i][j], t[i][j]);
			d[i][j] = max(d[i - 1][j - 1] + (nst_nt4_table[(int)text[i - 1]] == nst_nt4_table[(int)ref[j - 1]] ? opt->match : -opt->miss), max(r[i][j], t[i][j]));
		}
	}

	i = m;j = n;

	x = 0;
    {
        if (d[i][j] == r[i][j])
        {
            cigar[x].c = 'D';
            cigar[x].l = 1;
			j--;
		}
		else if (d[i][j] == t[i][j])
        {
            cigar[x].c = 'I';
            cigar[x].l = 1;
			i--;
		}
		else
        {
            if(d[i][j]>d[i-1][j-1])
            {
                cigar[x].c = 'M';
                cigar[x].l = 1;
            }
            else
            {

                cigar[x].c = 'X';
                cigar[x].l = 1;
            }
			i--, j--;
		}
    }

	while (i > 0 || j > 0)
    {
		if (d[i][j] == r[i][j])
        {
            if(cigar[x].c=='D') cigar[x].l++;
            else
            {
                x++;
                if(x>=MAX_CIGAR_BUF) {return -10000;}
                cigar[x].c = 'D';
                cigar[x].l = 1;
            }
			j--;
		}
		else if (d[i][j] == t[i][j])
        {
            if(cigar[x].c=='I') cigar[x].l++;
            else
            {
                x++;
                if(x>=MAX_CIGAR_BUF) {return -10000;}
                cigar[x].c = 'I';
                cigar[x].l = 1;
            }
			i--;
		}
		else
        {
            if(d[i][j]>d[i-1][j-1])
            {
                if(cigar[x].c=='M') cigar[x].l++;
                else
                {
                    x++;
                    if(x>=MAX_CIGAR_BUF) {return -10000;}
                    cigar[x].c = 'M';
                    cigar[x].l = 1;
                }
            }
            else
            {
                if(cigar[x].c=='X') cigar[x].l++;
                else
                {
                    x++;
                    if(x>=MAX_CIGAR_BUF) {return -10000;}
                    cigar[x].c = 'X';
                    cigar[x].l = 1;
                }
            }
			i--, j--;
		}
	}
	x++;
	if(x>=MAX_CIGAR_BUF) {return -10000;}

	char c;
    int l;
	for(i = 0;i<x/2;i++)
    {
        c = cigar[i].c;
        l = cigar[i].l;

        cigar[i].c = cigar[x-i-1].c;
        cigar[i].l = cigar[x-i-1].l;

        cigar[x-i-1].c = c;
        cigar[x-i-1].l = l;
    }

    (*cigar_num) = x;
    for(i = 0;i<x;i++)
    {
        if((cigar[i].c=='D')&&(cigar[i].l>opt->change_length)) cigar[i].c='N';
        if((cigar[i].c=='I')&&(cigar[i].l>20)) d[m][n]=-10000;
    }

	return d[m][n];
}
int Area_DP(char *chr,unsigned int start,unsigned int end,char* string,int length,struct cigar_t *cigar,int *cigar_num,struct Splice_DP_t *DP)
{
    int GAP = 4;
    int MaxPenalty = -65536;

    char *text = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
    char *ref = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
    struct DPC_t *cigarD = DP->cigarD;
    int DN = 0;

    int i, j,x;

    int n = end-start+1;
    int  m = length;

    if((n>=MAX_DP_AREA)||(m>=MAX_DP_AREA))
    {
        (*cigar_num) = 0;
        free(ref);
        free(text);
        return -10000;
    }

    for(i = 0;i<n;i++)
        ref[i] = chr[start+i];
    ref[n] = '\n';

    for(i = 0;i<m;i++)
        text[i] = string[i];
    text[m] = '\0';

	int** r =  DP->r;
	int** t =  DP->t;
	int** d =  DP->d;

	// initialization
	r[0][0] = t[0][0] = d[0][0] = 0;
	for (i = 1; i <= m; i++)
	{
		r[i][0] = MaxPenalty;
		d[i][0] = t[i][0] = -(GAP + (i-1)*opt->gap);
	}

	for (j = 1; j <= n; j++)
	{
		t[0][j] = MaxPenalty;
		d[0][j] = r[0][j] = -min(opt->splice,(GAP + (j-1)*opt->gap));//-(GAP + (j-1)*opt->gap);//-min(opt->splice,(GAP + (j-1)*opt->gap));
	}

	for (i = 1; i <= m; i++)
	{
		for (j = 1; j <= n; j++)
		{
			r[i][j] = max(r[i][j - 1] - opt->gap, d[i][j - 1] - GAP);
			t[i][j] = max(t[i - 1][j] - opt->gap, d[i - 1][j] - GAP);
			//if (s1[i - 1] == 'N' || s2[j - 1] == 'N') s[i][j] = max(s[i - 1][j - 1] + 1, r[i][j], t[i][j]);
			//else s[i][j] = max(s[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 1 : -1), r[i][j], t[i][j]);
			d[i][j] = max(d[i - 1][j - 1] + (nst_nt4_table[(int)text[i - 1]] == nst_nt4_table[(int)ref[j - 1]] ? opt->match : -opt->miss), max(r[i][j], t[i][j]));
		}
	}

	i = m;j = n;
	DN = 0;
	while (i > 0 || j > 0)
    {
		if (d[i][j] == r[i][j])
        {
            cigarD[DN].cigar = 'D';
            cigarD[DN].score = d[i][j];
            cigarD[DN].m = i;
            cigarD[DN].n = j+start-1;
            DN++;
			j--;
		}
		else if (d[i][j] == t[i][j])
        {
			cigarD[DN].cigar = 'I';
            cigarD[DN].score = d[i][j];
            cigarD[DN].m = i;
            cigarD[DN].n = j+start-1;
            DN++;
			i--;
		}
		else
        {
            if(d[i][j]>d[i-1][j-1]) cigarD[DN].cigar = 'M';
            else cigarD[DN].cigar = 'X';
            cigarD[DN].score = d[i][j];
            cigarD[DN].m = i;
            cigarD[DN].n = j+start-1;
            DN++;
			i--, j--;
		}
	}

    x = 0;
    cigar[x].c = cigarD[DN-1].cigar;
    cigar[x].l = 0;
    for(i=DN-1;i>=0;i--)
    {
        if(cigar[x].c==cigarD[i].cigar)
            cigar[x].l++;
        else
        {
            x++;
            cigar[x].c = cigarD[i].cigar;
            cigar[x].l = 1;
        }
    }
    if(cigar[x].l != 0)x++;
    (*cigar_num) = x;
    for(i = 0;i<x;i++)
    {
        if((cigar[i].c=='D')&&(cigar[i].l>opt->change_length)) cigar[i].c='N';
    }
    free(ref);
    free(text);
	return cigarD[0].score;
}
int Splice_DP(char *chr,unsigned int start,unsigned int end,char* string,int length,struct cigar_t *cigar,int *cigar_num,struct Splice_DP_t *DP)
{
    int GAP = 4;
    //int GAP2 = 24;
    int MaxPenalty = -65536;

    int GTAG = 2*4*opt->match;

    char *text = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
    char *ref = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
    struct DPC_t *cigarD = DP->cigarD;
    int DN = 0;
    //struct DPC_t *cigarA = DP->cigarA;
    //int AN = 0;

    int i, j;

    int area = end-start+1;
    int  m = length;
    int n = min(max(length*(1.1),length+10),area);

    if((n>=MAX_DP_AREA)||(m>=MAX_DP_AREA))
    {
        (*cigar_num) = 0;
        free(ref);
        free(text);
        return -10000;
    }

    for(i = 0;i<n;i++)
        ref[i] = chr[start+i];
    ref[n] = '\n';

    for(i = 0;i<m;i++)
        text[i] = string[i];
    text[m] = '\0';

	int** r = DP->r;
	int** t =  DP->t;
	int** d =  DP->d;
	int** a =  DP->a;

	struct max_t *maxD = DP->max_d;
	struct max_t *maxA = DP->max_a;
	// initialization
	maxD[0].m = maxD[0].score =r[0][0] = t[0][0] = d[0][0] = 0;
	for (i = 1; i <= m; i++)
	{
	    maxD[i].score = MaxPenalty;maxD[i].m = 0;
		r[i][0] = MaxPenalty;
		d[i][0] = t[i][0] = -(GAP + (i-1)*opt->gap);
	}

	for (j = 1; j <= n; j++)
	{
		t[0][j] = MaxPenalty;
		d[0][j] = r[0][j] = -min(opt->splice,(GAP + (j-1)*opt->gap));
	}

	for (i = 1; i <= m; i++)
	{
		for (j = 1; j <= n; j++)
		{
			r[i][j] = max(r[i][j - 1] - opt->gap, d[i][j - 1] - GAP);
			t[i][j] = max(t[i - 1][j] - opt->gap, d[i - 1][j] - GAP);
			d[i][j] = max(d[i - 1][j - 1] + (nst_nt4_table[(int)text[i - 1]] == nst_nt4_table[(int)ref[j - 1]] ? opt->match : -opt->miss), max(r[i][j], t[i][j]));
            if(d[i][j]>maxD[i].score) {maxD[i].score = d[i][j];maxD[i].m = j;}
		}
	}

	for(i = 0;i<n;i++)
        ref[i] = chr[end-i];
    ref[n] = '\n';

    for(i = 0;i<m;i++)
        text[i] = string[m-i-1];
    text[m] = '\0';

    r = DP->r1;
    t = DP->t1;

	maxA[0].m = maxA[0].score = r[0][0] = t[0][0] = a[0][0] = 0;
	for (i = 1; i <= m; i++)
	{
	    maxA[i].score = MaxPenalty;maxA[i].m = 0;
		r[i][0] = MaxPenalty;
		a[i][0] = t[i][0] = MaxPenalty;//-(GAP + (i-1)*opt->gap);
	}

	for (j = 1; j <= n; j++)
	{
		t[0][j] = MaxPenalty;
		a[0][j] = r[0][j] = -(GAP + (j-1)*opt->gap);
	}

	for (i = 1; i <= m; i++)
	{
		for (j = 1; j <= n; j++)
		{
			r[i][j] = max(r[i][j - 1] - opt->gap, a[i][j - 1] - GAP);
			t[i][j] = max(t[i - 1][j] - opt->gap, a[i - 1][j] - GAP);
			a[i][j] = max(a[i - 1][j - 1] + (nst_nt4_table[(int)text[i - 1]] == nst_nt4_table[(int)ref[j - 1]] ? opt->match : -opt->miss), max(r[i][j], t[i][j]));
            if(a[i][j]>maxA[i].score) {maxA[i].score = a[i][j];maxA[i].m = j;}
		}
	}

	int score = MaxPenalty;
	int max_x = 0,max_y = 0;

	for(i = 0;i<=m;i++)
    {
        if((maxD[i].score+maxA[m-i].score)>score) {score = (maxD[i].score+maxA[m-i].score);max_x = i;max_y = m-i;}
        else if((nst_nt4_table[(int)chr[(maxD[i].m+start-1+1)]]==1)&&(nst_nt4_table[(int)chr[(maxD[i].m+start-1+2)]]==3)
                &&(nst_nt4_table[(int)chr[(end-maxA[m-i].m+1-2)]]==0)&&(nst_nt4_table[(int)chr[(end-maxA[m-i].m+1-1)]]==1)
                &&((maxD[i].score+maxA[m-i].score+GTAG)>score))
            {score = (maxD[i].score+maxA[m-i].score+GTAG);max_x = i;max_y = m-i;}
        else if((nst_nt4_table[(int)chr[(maxD[i].m+start-1+1)]]==2)&&(nst_nt4_table[(int)chr[(maxD[i].m+start-1+2)]]==3)
                &&(nst_nt4_table[(int)chr[(end-maxA[m-i].m+1-2)]]==0)&&(nst_nt4_table[(int)chr[(end-maxA[m-i].m+1-1)]]==2)
                &&((maxD[i].score+maxA[m-i].score+GTAG)>score))
            {score = (maxD[i].score+maxA[m-i].score+GTAG);max_x = i;max_y = m-i;}
    }

    r = DP->r;
    t = DP->t;
    i = max_x;j = maxD[max_x].m;
    while (i > 0 || j > 0)
    {
		if (d[i][j] == r[i][j])
        {
            cigarD[DN].cigar = 'D';
            DN++;
			j--;
		}
		else if (d[i][j] == t[i][j])
        {
			cigarD[DN].cigar = 'I';
            DN++;
			i--;
		}
		else
        {
            if(d[i][j]>d[i-1][j-1]) cigarD[DN].cigar = 'M';
            else cigarD[DN].cigar = 'X';
            DN++;
			i--, j--;
		}
	}
	int x = 0;
    cigar[x].c = cigarD[DN-1].cigar;
    cigar[x].l = 0;
    for(i=DN-1;i>=0;i--)
	//for(i=max_x;i<DN;i++)
    {
        if(cigar[x].c==cigarD[i].cigar)
            cigar[x].l++;
        else
        {
            x++;
            cigar[x].c = cigarD[i].cigar;
            cigar[x].l = 1;
        }
    }
    if(cigar[x].l != 0)x++;

    cigar[x].c = 'N';
    if((max_x!=0)&&(max_y!=0))cigar[x].l = (end-maxA[max_y].m+1)-(maxD[max_x].m+start-1)-1;
    else if(max_x==0)cigar[x].l = (end-maxA[max_y].m+1)-start;
    else if(max_y==0)cigar[x].l = end-(maxD[max_x].m+start-1);

    if(cigar[i].l<8){cigar[i].c='D';}

    if(cigar[x].l>0) x++;
    else if(cigar[x].l<0) {free(ref);free(text);return Area_DP(chr,start,end,string,length,cigar,cigar_num,DP);}

    r = DP->r1;
    t = DP->t1;
    i = max_y;j = maxA[max_y].m;
    DN = 0;
	while (i > 0 || j > 0)
    {
		if (a[i][j] == r[i][j])
        {
			cigarD[DN].cigar = 'D';
            DN++;
			j--;
		}
		else if (a[i][j] == t[i][j])
        {
			cigarD[DN].cigar = 'I';
            DN++;
			i--;
		}
		else
        {
			if(a[i][j]>a[i-1][j-1]) cigarD[DN].cigar = 'M';
            else cigarD[DN].cigar = 'X';
            DN++;
			i--, j--;
		}
	}

    cigar[x].c = cigarD[0].cigar;
    cigar[x].l = 0;
    for(i = 0;i<DN;i++)
    {
        if(cigar[x].c==cigarD[i].cigar)
            cigar[x].l++;
        else
        {
            x++;
            cigar[x].c = cigarD[i].cigar;
            cigar[x].l = 1;
        }
    }
    if(cigar[x].l != 0)x++;
    (*cigar_num) = x;

    score = 0;
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

    free(ref);
    free(text);
	return score;
}
int tail_seed(int seed_length,int chr,char strand,uint64_t pos,char *seq,int length,int mode,struct seed_t *seed_o)//0 front 1 back
{
    //if(length>200){(*cigar_num) = 0;return 1;}
    int i,j,k;
    //struct seed_t seed_r;
    struct seed_t seed[101];
    int seed_num = 0;

    unsigned int chr_start = opt->chr->list[chr].start_site;

    if(mode==0)
    {
        int l = min(10,length-seed_length);
        for(j = 0-opt->change_length;j<=opt->change_length;j++)//for(j = -8;j<=8;j++)//
        {
            for(i = 0;i<l;i++)
            {
            if(strand==1)
            {
                k = 0;
                while((i+k<length)&&(3-nst_nt4_table[(int)seq[length-i-k-1]]==nst_nt4_table[(int)opt->chr->list[chr].seq[(opt->l_pac<<1)-(pos-j-i-k-1)-1-chr_start]]))
                    k++;
            }
            else
            {
                k = 0;
                while((i+k<length)&&(nst_nt4_table[(int)seq[length-i-k-1]]==nst_nt4_table[(int)opt->chr->list[chr].seq[pos-j-i-k-1-chr_start]]))
                    k++;
            }

            if(k>=seed_length)
            {
                seed[seed_num].pos = pos-j-k-i;
                seed[seed_num].start = length-i-k;
                seed[seed_num].length = k;
                seed[seed_num].last = -1;
                seed[seed_num].abs = seed[seed_num].pos-seed[seed_num].start;

                if(seed[seed_num].pos+seed[seed_num].length>pos)
                {
                    seed[seed_num].length-=seed[seed_num].pos+seed[seed_num].length-pos;
                }
                if(seed[seed_num].length>=seed_length)
                {
                    seed[seed_num].score = seed[seed_num].length;
                    seed_num++;
                    //if(seed_r.abs!=pos-1) seed_r.score-=opt->gap;
                    //if(seed_r.abs>pos-1) seed_r.score-=seed_r.abs-pos+1;
                    //insert_seed(seed,&seed_num,100,seed_r);

                }
            }
            i+=k;
            }
        }
    }
    else if(mode==1)
    {
        int l = min(10,length-seed_length);
        for(j = 0-opt->change_length;j<=opt->change_length;j++)//for(j = -8;j<=8;j++)//
        {
            for(i = 0;i<l;i++)
            {
                if(strand==1)
                {
                    k = 0;
                    while((i+k<length)&&(3-nst_nt4_table[(int)seq[i+k]]==nst_nt4_table[(int)opt->chr->list[chr].seq[(opt->l_pac<<1)-(pos+1+j+i+k)-1-chr_start]]))
                        k++;
                }
                else
                {
                    k = 0;
                    while((i+k<length)&&(nst_nt4_table[(int)seq[i+k]]==nst_nt4_table[(int)opt->chr->list[chr].seq[pos+1+j+i+k-chr_start]]))
                        k++;

                }
                if(k>=seed_length)
                {
                    seed[seed_num].pos = pos+1+j+i;
                    seed[seed_num].start = i;
                    seed[seed_num].length = k;
                    seed[seed_num].last = -1;
                    seed[seed_num].abs = seed[seed_num].pos-seed[seed_num].start;

                    if(seed[seed_num].pos<=pos)
                    {
                        seed[seed_num].length-=pos-seed[seed_num].pos+1;
                        seed[seed_num].start+=pos-seed[seed_num].pos+1;
                        seed[seed_num].pos+=pos-seed[seed_num].pos+1;
                    }
                    if(seed[seed_num].length>=seed_length)
                    {
                        seed[seed_num].score = seed[seed_num].length;
                        seed_num++;
                        //if(seed_r.abs!=pos+1) seed_r.score-=opt->gap;
                        //if(seed_r.abs<pos+1) seed_r.score-=pos+1-seed_r.abs;
                        //insert_seed(seed,&seed_num,100,seed_r);
                    }
                }
                i+=k;
            }
        }
    }
    int max_o,max_s;

    for(i = 0; i<seed_num; i++)
    {
        seed[i].score = seed[i].length*(opt->match);

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

                    if(seed[i].abs<seed[j].abs)
                        max_s -=(seed[j].abs-seed[i].abs)*opt->match+4+(seed[j].abs-seed[i].abs-1)*opt->gap;
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
    }
        //回溯
    max_o = -1;
    max_s = -1000;

    for(i = seed_num-1; i>=0; i--)
    {
        if(seed[i].score>=max_s)
        {
            max_o = i;
            max_s = seed[i].score;
        }
    }
    //回溯
    if(max_o==-1)
    {
        return 0;
    }
    else
    {
        int x = 0;
        while(max_o!=-1)
        {
            seed_o->pos = seed[max_o].pos;
            seed_o->start = seed[max_o].start;
            seed_o->length = seed[max_o].length;
            seed_o->abs = 0;
            seed_o->last = -1;
            seed_o->score = 0;
            x++;
            seed_o++;

            if(seed[max_o].last!=-1)
                max_o = seed[max_o].last;
            else max_o = -1;
        }
        //if(x>1)
            //printf("1\n");
        return x;
    }
}


